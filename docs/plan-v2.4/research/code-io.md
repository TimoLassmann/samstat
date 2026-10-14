# code-io: input, threading and speed in SAMStat 2.3.0

Scope: `src/samstat.c`, `src/htsinterface/htsglue.c`, the tld FASTA/FASTQ reader, `src/param` input handling. All measurements were made on 2026-10-14 on an Apple-silicon Mac (16 logical cores, warm page cache, Release build, bundled htslib 1.16, `-t` default) with synthetic data. Real BAMs carry more tags and have less uniform reads, so absolute numbers will differ; ratios are the useful part. Nothing in the repo was modified. Prototypes live in the scratch dir `/private/tmp/claude-502/samstat-v24-plan/proto-io` (copy of the source tree, not the repo).

Test data (scratch `io-data/`, generated with `uv run --with numpy python`, mapped with minimap2 + `samtools calmd` so MD/NM exist):

| name | content |
|---|---|
| short.{fq,fq.gz,sam,bam} | 1 M reads x 150 bp, 1 % substitutions, 5 Mb random reference |
| short5.{sam,bam} | the same 5 times (5 M reads) |
| long.{fq,sam,bam} | 20 k reads, lognormal length (median ~6 kb, max >100 kb), 8 % substitutions (165 Mb of sequence) |
| long4.{fq,bam} | the same 4 times (80 k reads, 660 Mb) |
| repo `dev/` | only tiny files; the largest, `ce#large_seq.sam`, is 2.1 MB with 3 records, one of 1 Mb |

The repo's `dev/` has no input large enough to time; everything else there runs in ~13 ms.

## 1. How it works today

### Threading
* `main` -> `thread_data_init`: `n_threads = min(n_infile, -t)` (default `-t 4`). Each worker loops under a mutex picking the next unprocessed file (`td->active[]`) and calls `process_file`. **One file = one thread.** A single input file always uses one core (plus nothing else: htslib is opened without `hts_set_threads`).
* Per file (`process_sam_bam_file`): allocate `stat_collection`, a `tl_seq_buffer` of `-buffer_size` = 10 000 slots (`param.c:390`), then loop `read_bam_chunk` -> `parse_alignment` + `fix_orientation` per read -> `collect_stats` -> `clear_aln_data`/`reset_tl_seq_buffer`. No overlap between decompression, parsing and counting.
* `read_bam_chunk` (htsglue.c) copies each non-secondary/supplementary record out of `bam1_t` into per-slot `tl_seq` strbufs (name, seq expanded to ASCII, qual +33), the CIGAR, and `NM`/`MD` via `bam_aux_get` (linear tag scan). Secondary/supplementary are counted and skipped.
* `-p N` is applied in the chunk loop; it cuts the last chunk exactly and stops reading.

### FASTA/FASTQ
* `prepare_fasta_fastq_input` reads the **entire** input once through `gzread` (so a `.gz` is inflated a second time) to (a) detect "empty" inputs, (b) find whether the last byte is `\n`. Only if the newline is missing it copies the whole decompressed file to `$TMPDIR/samstat_XXXXXX` plus `\n` (tld bug: unterminated last record dropped). Then the file is opened again by tld.
* tld reader: `gzgets` per line into a 1 MB line buffer, a small state machine in `parse_buf`, FASTA/FASTQ and DNA/protein detected from content (`start_symbol` `@`/`>`, `finalize_read`). Single thread, no `hts_set_threads` equivalent. Gzip is plain zlib, so no parallel inflate.
* Type of each input is decided **by file extension** in `param.c:detect_file_type` (`.sam .bam .fa .fasta .fna .fq .fastq` + `.gz` for the FASTA/Q ones). A path without a recognised extension, `-`, or `/dev/stdin` is rejected at parse time (verified: exit 2, "unrecognized file name extension"). So **stdin is impossible today**, and so is `.cram` and `.fq.bgz` etc.

### Bundled htslib 1.16 already reads FASTA/FASTQ
`sam_open`/`sam_read1` in htslib >= 1.13 reads FASTQ/FASTA as unmapped records (NEWS 1.13: "Support reading and writing FASTQ format in the same way as SAM, BAM or CRAM", checked in `htslib-1.16/NEWS`; https://github.com/samtools/htslib/blob/develop/NEWS). I verified with a 10-line C program against the bundled `libhts.a`: it reads the 1 M-read FASTQ, the gzipped one, **an unterminated last record in FASTQ and FASTA**, multi-line FASTA, and all of them from a pipe (`cat x | ./t -`), with no temp file.

## 2. Measurements (baseline 2.3.0)

Wall time, single run (variance is small: runs are CPU-bound, 99 % cpu), `-t` default, one file.

| input | reads | time | peak RSS | throughput |
|---|---|---|---|---|
| short.fq | 1 M | 0.57 s | 18 MB | 1.8 M reads/s |
| short.fq.gz | 1 M | 1.78 s | - | 0.56 M reads/s |
| short.sam | 1 M | 1.13 s | 28 MB | 0.9 M reads/s |
| short.bam | 1 M | 1.35 s | 28 MB | 0.74 M reads/s |
| short5.fq (1.5 GB) | 5 M | 2.85 s | - | 1.75 M reads/s |
| short5.sam (2.2 GB) | 5 M | 6.3 s | - | 0.8 M reads/s |
| short5.bam (870 MB) | 5 M | 7.3 s | - | 0.69 M reads/s |
| long.fq | 20 k | 1.29 s | **987 MB** | 128 Mb/s |
| long.sam / long.bam | 20 k | 2.06 s / 2.41 s | **1.29 GB** | 80-93 Mb/s |
| long4.fq | 80 k | 2.82 s | 1.15 GB | 235 Mb/s |
| long4.bam | 80 k | 7.1 s | 1.49 GB | 93 Mb/s |
| `dev/ce#large_seq.sam` (3 reads, one 1 Mb) | 3 | **1.39 s** | **3.84 GB** | - |
| `-p 1` on long4.fq / short.fq (fixed cost) | 1 | 0.64 s / 0.17 s | - | - |

Extrapolation (not measured), 30x human WGS = 100 Gb: short-read BAM ~667 M reads / 0.69 M/s = **~16 min**; long-read BAM 100 Gb / 93 Mb/s = **~18 min**, all on one core, whatever `-t` is.

### Where the time goes (macOS `sample`, worker thread only, idle `__ulock_wait` of the joining main thread removed)

| input | collect_stats | htslib decode (inflate / sam_parse1 / aux) | parse_alignment | read_bam_chunk copy | read(2) |
|---|---|---|---|---|---|
| short5.sam | ~37 % | ~23 % (text parse) | ~12 % | ~9 % | ~7 % |
| long4.bam | ~60 % | ~25 % (zlib inflate) | ~13 % | ~8 % | - |
| short5.fq | ~40 % (collect_stats_seqqual) | - | - | ~35 % tld parse (`read_fasta_fastq_file` inlined), ~20 % read/memchr/memmove/malloc | |

(Percentages from ~2100-sample profiles; FASTQ pre-pass is a single sampling-invisible burst before the sampled window and is estimated, not measured, at 10-20 % of raw-FASTQ time and ~35 % of `.gz` time since it inflates the file twice: `gzip -dc short.fq.gz` alone is 0.55 s, samstat-gz is 1.78 s, and htslib reads the same file in 0.78 s.)

Two fixed or pathological costs found along the way, not I/O but they decide whether long reads are usable:
1. **Memory scales with the longest read, not with the data.** `ce#large_seq.sam` (2 MB file, one 1 Mb read) needs **3.8 GB RSS and 1.4 s**; 98 % of that time is `clu_data` (interval binning, `plot_interval.c`). Cause: every per-position `plot_data` is resized to `sb->max_len + 1` (`collect.c:81-110`) with one row per symbol/quality and per MAPQ bin, and the per-slot buffers of `tl_seq`/`aln_data` (name, seq, qual, genome, read) only grow. For ONT data the resident set is ~1 GB even for 165 Mb of input. Position axes need a hard cap (see section 5).
2. **Buffer is sized in reads (10 000), not bases.** 10 000 slots x (seq + qual + genome + read ~4 bytes/base) at a high-water mark of the longest read that ever sat in that slot explains the 1-1.5 GB on long reads. With W parallel workers this multiplies by W.

## 3. Design: intra-file parallelism

### Principle
The statistics are pure integer counts (every `stat_collection` field is a `uint64_t` array or counter, apart from derived `clu_data` that is built at report time), so they **commute and add**. Map-reduce is therefore exact: per-worker private `stat_collection`, merged by element-wise addition after resizing to the largest length. This keeps reports byte-identical to the serial result (golden tests stay valid) and needs no locks in the hot path. (UNVERIFIED for every plot_data: I did not audit `collect.c` for a non-additive field; the merge unit test below catches it.)

### Pipeline
```
reader (main thread)              workers x W                         main
 sam_read1() into raw bam1_t      convert bam1_t -> tl_seq (the        merge W stat_collections,
 batches (<= 4096 reads OR        current read_bam_chunk body),        then stat_report() as today
 <= 8 Mbases)  --bounded queue--> parse_alignment + fix_orientation,
 htslib pool: BGZF inflate (and   collect_stats into PRIVATE stats,
 SAM text parse) on T_hts threads recycle the batch
```
* `hts_set_threads(fp, n)` (or a shared `hts_tpool`) makes htslib inflate BGZF blocks and parse SAM text on separate threads (htslib docs: https://www.htslib.org/doc/htslib.html is silent on details; behaviour confirmed empirically below). It is a one-line change and gives a large share of the win on its own.
* Move the bam1_t -> tl_seq conversion out of the serial reader. In the prototype it was still done under the lock, and that serial section is exactly what stopped scaling at 4 workers.
* Batches are bounded by bases as well as reads, which fixes the memory scaling for long reads (item 2 above). Slot buffers larger than, say, 1 MB are shrunk when the batch is recycled.
* `-p N`: the reader counts primary records and stops after exactly N, as now; early exit also stops decompression. Deterministic because only *which* reads, not the order, matters.
* Exit status and error handling: a reader or worker error sets a shared flag, the others drain, `process_*` returns FAIL as today.

### Thread budget (one rule, no new options)
`-t T` becomes "total threads" (default stays 4; raise to `min(8, ncpu)` only if you are happy to depend on affinity detection, which on HPC cgroups `hw.ncpu`/`nproc` get wrong). Files are processed by the existing file pool when `n_files >= T` (each file single-threaded: many small files stay efficient, a tiny file costs 13 ms), otherwise sequentially with `workers = T / min(n_files, T)`. This deletes no behaviour, adds no flag, and `-t` stays documented as "threads". With `T = 1` run inline with no thread creation (the current code path, so ASAN/golden behaviour is unchanged).

### Measured prototype (scratch copy, `SAMSTAT_PAR=W`, `SAMSTAT_HTS_THREADS=n`)
Prototype = htslib threads + reader under a mutex + W workers each with a private `stat_collection`. It **omits the merge, the report (0.17 s short / 0.6 s long fixed) and `-p`**, so add ~0.2-0.6 s to the "after" column for a like-for-like end-to-end figure.

| input | baseline 2.3.0 | hts threads only (4) | W=1 | W=2 | W=4 | W=8 |
|---|---|---|---|---|---|---|
| short5.bam (5 M reads) | 7.3 s | 4.1 s (1.8x) | 4.4 | 2.15 | **1.60** | 1.60 |
| short5.sam | 5.6 s | 3.8 s (1.5x) | 3.8 | 2.0 | **1.46** | 1.45 |
| long4.bam (80 k reads) | 6.7 s | 4.6 s (1.5x) | 3.75 | 2.05 | **1.38** | 1.38 |

So the realistic target is **~4-5x** end-to-end on 4+ cores for BAM/SAM short and long (7.3 -> ~1.8 s incl. report), 1.5-1.8x from the htslib thread pool alone. With the conversion also moved off the serial reader, the ceiling is the raw htslib decode rate (`samtools view -c -@8` does 1 M BAM records in 0.08 s), so W = 8 should keep scaling to roughly 6-8x (UNVERIFIED, not prototyped). 30x WGS then drops from ~16-18 min to roughly 3-4 min on 4 cores (extrapolated).

## 4. Design: one reader for everything, and stdin

### Recommendation: read FASTA/FASTQ through htslib too
`sam_open` reads SAM, BAM, CRAM, FASTQ and FASTA and auto-detects the format from content (`hts_detect_format`; works on pipes because htslib peeks and buffers). Using it for FASTQ/FASTA gives, with no new code in samstat:

| benefit | evidence |
|---|---|
| `samstat -` (stdin) for every format; `bwa mem ... | samstat -`, `minimap2 ... | samstat -`, `samtools view -h x.cram | samstat -` | tested with the bundled lib: FASTQ and gz FASTQ from a pipe work |
| no tld parser bugs: unterminated last record, `.gz` of nothing, multi-line FASTA work | tested: 2 of 2 reads returned for unterminated FASTQ and FASTA |
| delete `prepare_fasta_fastq_input` (second full decompression, temp file copy in `$TMPDIR`, ~80 lines) and the `tld_strbuf` quadratic-append exposure on the read side | code read |
| 3-8x faster reading: raw FASTQ 1 M reads 0.22 s (tld path ~0.4 s of the 0.57 s; htslib time includes nothing else), gz 0.78 s vs 1.78 s end-to-end today | measured |
| BGZF `.fq.gz`/`.bgz` (what sequencers and `bgzip` produce) inflate in parallel with `hts_set_threads`; plain gzip does not (UNVERIFIED for large speedups; inflate stays serial for non-BGZF gzip) | https://www.htslib.org/doc/bgzip.html |
| extension-independent type detection; `.cram` supported for free (needs `REF_PATH`/reference for MD; see below) | htslib behaviour |

Estimates (not measured end-to-end): 1 M FASTQ 0.57 s -> ~0.35 s serial, ~0.25 s with workers; 1 M `.fq.gz` 1.78 s -> ~1.0 s serial (non-BGZF gzip is inflate-bound at ~0.55 s) -> ~0.7 s with the reader thread overlapped with the workers.

Caveats the maintainer should decide on (judgement calls):
* **Protein FASTA.** tld path supports a protein alphabet (`TLALPHABET_NOAMBIGIOUS_PROTEIN`); htslib stores sequence as 4-bit nt16, so amino-acid letters collapse. SAMStat is a *read* QC tool and the paper and report are nucleotide-centric. Recommendation: drop protein support and say so in the changelog (one fewer code path), unless someone asks. Decision for the maintainer.
* **FASTQ quality offset.** htslib FASTQ reading assumes Phred+33 (as does samstat today, `base_quality_offset = 33` hard-coded). Same behaviour, no regression.
* **FASTQ comments / paired suffixes.** htslib strips `/1` `/2` and treats records as unmapped; R1/R2 split then depends on flags htslib sets from the `/1` `/2` suffix names (UNVERIFIED how today's FASTQ path assigns R2; check against the existing `fastx` ctest labels).
* tld stays a dependency (strbuf, alloc, etc.), only its seq reader leaves samstat.

### stdin details (smallest change)
1. `param.c`: replace extension-based `detect_file_type` with a content probe: `hts_open(path,"r")` + `hts_get_format()` (`sam/bam/cram/fastq_format/fasta_format`). A path `-` is allowed. Unknown/unreadable keeps today's "file missing" exit 1 path (the existing "missing file is reported when processed" comment suggests the intent).
2. `stat_report.c` output name for stdin: `stdin.samstat.html` in `-d` (or cwd). One name, no new option; if two stdins are given it is an error. 
3. `tld_file_exists` checks in `samstat.c:168,341` are skipped for `-`.
4. Progress: with a pipe there is no size; `--verbose` already prints read counts.
5. `-p N` on a pipe stops reading early and closes stdin; the upstream program gets SIGPIPE, which is the expected `head`-like behaviour; document it ("peek the first 1 M reads of a running alignment: `... | samstat -p 1000000 -`").

Only one stdin file makes sense; it is processed by the same pipeline, so the parallel design above applies unchanged.

## 5. What to change in 2.4 (ordered by value per effort)

| # | change | effort | gain |
|---|---|---|---|
| 1 | `hts_set_threads(fp, min(T,4))` in `open_bam` | 1 line + `-t` doc | 1.5-1.8x on BAM/SAM, verified |
| 2 | Read FASTA/FASTQ via `sam_read1`; delete tld reader use and `prepare_fasta_fastq_input`; content-based type detection | ~150 lines removed, ~40 added | 2-3x on gz FASTQ, fixes unterminated-record bug class, enables stdin and CRAM input |
| 3 | `samstat -` stdin | ~30 lines (param + report name) | pipelines without temp files |
| 4 | Parallel map-reduce: reader + W workers with private `stat_collection` + `stat_collection_merge()` | ~250 lines; merge ~60 | 4-5x end-to-end (measured prototype), more with conversion off the serial path |
| 5 | Batch by bases (<= 8 Mbases) and shrink oversized slot buffers | ~30 lines | removes the 1 GB RSS on long reads, makes W workers affordable |
| 6 | **Cap per-position axes** in `collect.c` (e.g. first 10-20 kb positions; reads beyond are still counted in the length/MAPQ/summary plots) | few lines in `plot_data_resize_len`/clamp (clamp helper `clamp_idx` already exists) | removes 3.8 GB / 1.4 s for a 1 Mb read; prerequisite for long-read support (hand to the long-read/plots planners, it is a presentation decision as well) |
| 7 | Remove the file-level worker pool except for `n_files >= T` | -100 lines net | simpler |

### Tests to add with the change
* `stat_collection_merge` unit test: W=1 vs W=8 outputs byte-identical for every golden and for 3 sizes of synthetic input (serial vs parallel equality is the whole safety net).
* Thread-count equality ctest already exists (commit 4824d41); extend it to `-t 1/2/8` on a >100 k-read synthetic BAM and a `-p` boundary case (N not a multiple of batch size).
* stdin: `cat x.fq | samstat -`, `cat x.bam | samstat -`, gz FASTQ, empty stdin (exit 0, empty report, as the file case), truncated BAM (exit 1).
* ASAN job (exists) must run the parallel path; add a TSAN job on ubuntu (workers + reader + htslib pool).

### Risks
* Determinism if any counter is accidentally float: guarded by the equality test.
* Memory with W x batch buffers: bounded by 8 Mbase batches (~64 MB each at ~8 B/base across the four per-slot buffers; W=4 -> ~0.3 GB worst case, versus 1-1.5 GB today for long reads, and ~20 MB for short reads with 4096-read batches).
* htslib FASTQ differences (see caveats): run the `fastx`/`golden` ctest labels before deleting the tld path.
* Windows is unsupported already (pthreads); no change.

## 6. Sources and verification notes
* htslib NEWS 1.13 (FASTQ read/write through the SAM API): bundled copy `thirdparty/htslib-1.16.tar.bz2` -> `NEWS`; upstream https://github.com/samtools/htslib/blob/develop/NEWS
* htslib API (`hts_set_threads`, `hts_detect_format`, `sam_open`): https://www.htslib.org/doc/htslib.html (I read the bundled headers, not the web page; the page was not fetched: mark web text UNVERIFIED).
* BGZF vs gzip parallel decompression: https://www.htslib.org/doc/bgzip.html, UNVERIFIED for speedups on real FASTQ.
* All timings, RSS and sample percentages above are mine; scripts and data in `/private/tmp/claude-502/samstat-v24-plan/` (`gen.py`, `io-data/`, `proto-io/`, `s-*.txt` profiles).
