# perf-htslib: fast BAM/CRAM reading for SAMStat 2.4

Date: 2026-10-14. Scope: htslib threading, decode-vs-process cost, parallel designs, CRAM references, stdin, post-1.16 versions.

## 1. What SAMStat does today (src/htsinterface/htsglue.c)

- `sam_open(file,"r")`, `sam_hdr_read`, then a serial `sam_read1()` loop. No `hts_set_threads`, no thread pool, no `-T`/reference handling, no stdin path.
- Every non-secondary/supplementary record is copied into a `tl_seq` (name append, per-base `seq_nt16_str` loop, quality copy) and later re-parsed (CIGAR + MD) in `parse_alignment`.
- Parallelism is one pthread per input file (`samstat.c`). One large BAM = one core, and that core also decompresses BGZF.

Consequence: for one BAM the wall time is roughly (inflate) + (record decode) + (copy to tl_seq) + (CIGAR/MD parse) + (collect), all serial.

## 2. Where the time goes (typical, not measured here)

Rough, widely observed orders of magnitude; treat as UNVERIFIED for SAMStat's own build until benchmarked (see section 8):
- BGZF inflate (zlib) is single-thread ~150-300 MB/s of uncompressed output; libdeflate is ~1.5-2x faster. BAM at ~4:1 compression means inflate is usually the largest single cost of a bare `samtools view -c`-style loop.
- `bam_read1` decode into `bam1_t` is cheap (memcpy-like) compared with inflate.
- Per-record work in SAMStat (string appends, per-base char conversion, MD parsing, per-position counters) is likely comparable to or larger than inflate for short reads, and dominates for long reads. The base-by-base `seq_nt16_str` loop and the extra `tl_seq` copy are avoidable: parse straight from `bam1_t` (`bam_get_seq`, `bam_get_cigar`, `bam_aux_get(b,"MD")`) and accumulate.
- Documented throughput anchor: htslib's own docs/blog posts and samtools threading discussions report near-linear BGZF decompression scaling to ~4-8 threads, then flat once the single reader/consumer thread saturates. No authoritative MB/s table found; do not quote a number without a local benchmark. UNVERIFIED.

## 3. Cheapest win: BGZF/CRAM decompression threads

API (htslib `hts.h`): `hts_set_threads(htsFile*, int n)` or share a pool across files with `hts_tpool_init(n)` + `hts_set_thread_pool(fp, &(htsThreadPool){pool,0})`. Docs: https://www.htslib.org/doc/htslib.html , samtools `-@` in https://www.htslib.org/doc/samtools.html .
- One call after `sam_open` gives parallel inflate (read-ahead decoding blocks) for BAM, and parallel container decode for CRAM, while the caller keeps its serial `sam_read1` loop. Zero design change, zero merge logic.
- Expected benefit: removes inflate from the critical path; speedup bounded by the remaining serial per-record work (Amdahl). For a typical short-read BAM expect ~1.3-2x; more with libdeflate. UNVERIFIED, benchmark.
- Shared pool: with N input files and `-t T` workers, create ONE `hts_tpool` and attach it to every file, instead of each file spawning its own threads (avoids oversubscription). Since 1.22 `hts_tpool_worker_id()` exists (see NEWS link below).
- Gotcha: since 1.23.1 setting a thread pool on an htsFile a second time is disallowed. Set once, right after open.
- Gotcha: threaded `sam_read1` hang after EOF was fixed in 1.22; the bundled 1.16 may be affected. UNVERIFIED which versions exactly.
- Gotcha: the multi-thread CRAM iterator dropped long overlapping alignments (fixed 1.23). Irrelevant to sequential whole-file reads, relevant if regions are ever used.

## 4. Multi-threaded processing design (when decode threads are not enough)

Recommended minimal design (keeps SAMStat simple, deterministic output):
1. One reader thread: `sam_read1` loop (with htslib decode threads on, section 3) fills batches of N records (e.g. 4k-16k reads; for long reads batch by bases, ~32-64 Mb) into a ring/queue of `bam1_t` arrays (`bam_copy1` or swap pointers).
2. W workers pop batches, run parse_alignment + update a PER-WORKER accumulator (the existing collect structs, zero shared writes, no locks in the hot path).
3. At end, merge accumulators by elementwise sum (all SAMStat statistics are counts/histograms, so merge is trivially associative and the result is independent of batch order -> output identical for any `-t`). This matches the existing "thread-count equality" test intent.
4. Per-file parallelism (today) stays as the outer level; use a global budget: threads_total = -t, split between files and per-file workers, with one shared htslib pool.

Cost/complexity: ~200-300 lines. Order-dependence trap: `-p N` read sampling and R1/R2 pairing must be decided on the reader thread (deterministic count of records seen), not in workers.
Alternative that needs no reader/worker code at all: just `hts_set_threads(fp, T)` (section 3). Recommend doing this first, measuring, and only adding the worker pool if the profile shows per-record work > ~50% of wall time.
Region-split parallelism over an index (one thread per chromosome via `sam_itr_queryi`) is an option for indexed BAM but needs .bai/.crai, fails for unmapped reads and long-read uniformity, and complicates sampling. Not recommended.

## 5. CRAM reference handling

- CRAM stores bases relative to a reference; decoding needs it unless the file was written with `no_ref`/embedded reference. Lookup order in htslib: `@SQ M5`/`UR` via the REF_PATH and REF_CACHE environment variables (default REF_PATH historically fetched from EBI ENA), or an explicit FASTA via `hts_set_fai_filename(fp, "ref.fa")` (samtools `-T`/`--reference`). Docs: https://www.htslib.org/doc/reference-sequences.html and https://www.htslib.org/doc/htslib-s3-plugin.html (remote access context).
- Since 1.22 htslib NO LONGER fetches references from EBI by default, and ships a `ref-cache` tool (https://github.com/samtools/htslib/blob/develop/NEWS). So with a bundled 1.16 SAMStat would silently try network access; with >=1.22 it fails fast. Either way SAMStat should expose `--reference FILE` (maps to `hts_set_fai_filename`) and print a clear error if CRAM decode fails for lack of reference, naming REF_PATH and the flag.
- Current code has `.cram` in no extension handling in MultiQC (dropped per commit cb33670), and no reference option. If CRAM stays supported, one option + one error message is the whole requirement.
- MD tag: CRAM does not store MD; htslib regenerates MD/NM when the `CRAM_OPT_DECODE_MD` option is set (default behaviour differs by version/reference availability, UNVERIFIED). SAMStat's mismatch profile depends on MD, so for CRAM it must either set `hts_set_opt(fp, CRAM_OPT_DECODE_MD, 1)` or compute mismatches itself from the reference. Also true for BAMs lacking MD (many long-read aligners such as minimap2 omit MD unless `--MD`): the long-read plan needs a fallback (`cs`/`=`/`X` CIGAR ops, or reference-based mismatch calling). Flag to the planning agents handling long reads.
- CRAM decode is heavier per record than BAM (rANS/ transforms, reference diff) so threads matter more there; CRAM 3.1 codecs got auto-vectorisation in 1.19 and decoding speed-ups in 1.18; 1.24 moves more cram2bam work into worker threads.

## 6. stdin ("-") streaming

- `sam_open("-","r")` reads stdin; htslib auto-detects SAM/BAM/CRAM from magic bytes. BAM/BGZF from a pipe works fine; threads work too (inflate is stream-based). Since 1.22 bgzip's pipe block size was raised to avoid slowdowns reading from pipes; since 1.19 `hts_close` no longer really closes stdin/stdout. Source: https://github.com/samtools/htslib/blob/develop/NEWS
- CRAM from a pipe works for sequential reading when the reference is supplied (no seeking needed for a full pass). UNVERIFIED for all CRAM edge cases.
- Real blockers are in SAMStat, not htslib: (a) report/JSON naming uses the filename, so "-" needs a `--name` or default "stdin"; (b) the FASTA/FASTQ path copies input to a temp file to work around the tld unterminated-record bug, which a stream would need to keep (spool to temp, or fix tld upstream); (c) two-pass anything (e.g. -p sampling requiring file size) must be streaming-compatible.
- Use case value: `minimap2 ... | samtools view -b - | samstat -` and `samtools view -h x.cram | samstat -`. Reasonable cheap feature (BAM/SAM only first), low maintenance.

## 7. htslib versions after 1.16

Source for all: https://github.com/samtools/htslib/blob/develop/NEWS (fetched 2026-10-14; the fetched excerpt lacked 1.17 entries; 1.17 itself is UNVERIFIED). Latest release in that file: 1.24 (9 Jul 2026).
- 1.18 (Jul 2023): CRAM decoding speed-ups; faster whole-chromosome index lookups.
- 1.19 / 1.19.1 (Dec 2023 / Jan 2024): CRAM 3.1 codec auto-vectorisation; stdin/stdout no longer really closed by hts_close; 1.19.1 fixes aux data stored uncompressed.
- 1.20 (Apr 2024): faster B-array tag parsing; CRAM security fixes (heap overflow in RG:Z).
- 1.21 (Sep 2024): `bgzf_read()` faster, 10-15% faster decoding via higher parallelism; kputll speed-up (~5% BAM->SAM); nibble2base speed-up. Directly relevant: SAMStat reads via bgzf, and nibble2base speeds `bam_get_seq` expansion paths.
- 1.22 / 1.22.1 (May/Jul 2025): default CRAM output 3.1; no EBI reference fetch by default; `ref-cache` tool; threaded sam_read1 EOF hang fixed; `hts_tpool_worker_id()`; libdeflate-based CRC32 when available.
- 1.23 / 1.23.1 (Dec 2025 / Mar 2026): faster SAM header parse (many @SQ), mod-tag (MM/ML) parsing speed-up, multi-thread CRAM iterator fix; CRAM decoder hardening CVEs; second `hts_set_thread_pool` on one file disallowed.
- 1.24 (Jul 2026): parallel cram2bam moves more work to worker threads; CRAM v4.0 experimental support removed; synced-reader error checking.

Upgrade recommendation: bump the bundled htslib 1.16 -> 1.24 (or at minimum >=1.21/1.23.1). Gains: faster bgzf_read, CRAM speed and security fixes (CRAM CVEs in 1.20/1.23.1 matter because SAMStat parses untrusted user files), pipe fixes, no surprise network fetches. Risk: low (stable API; `bam_hdr_t` alias still valid). Needs: refresh thirdparty tarball and checksum, rerun the 314 tests and goldens (expect no golden change). Also check the bundled build enables libdeflate (configure `--with-libdeflate`; needs the library) since inflate dominates; in the ExternalProject build this is a cheap flag if libdeflate is available, otherwise zlib-ng/zlib. UNVERIFIED whether the current ExternalProject enables it.

## 8. Recommendations (ordered by value/effort)

1. Add `hts_set_threads` (shared `hts_tpool`) tied to `-t`. ~10 lines. Biggest single gain per line of code. Default `-t` should become min(4, cores) so users need no options.
2. Bump htslib to 1.24, build with libdeflate if present.
3. Parse straight from `bam1_t`; drop the `tl_seq` copy and per-base char loop for the BAM path (reduces both CPU and memory; long reads currently copy 3x).
4. Benchmark before building the reader/worker design: run `perf`/`time` on a 5-10 GB BAM (short-read and ONT) with -t 1/4/8; add worker pool only if per-record work remains >50% of wall time. Per-worker accumulators + sum merge keeps results thread-count-invariant.
5. `--reference FILE` for CRAM plus an explicit error message; set `CRAM_OPT_DECODE_MD` or compute mismatches without MD (also needed for long-read BAMs without MD).
6. `samstat -` for stdin (BAM/SAM), named "stdin" unless `--name`.

Explicitly not recommended: region/index-based parallelism, custom BGZF handling, mmap tricks, GPU.

## 9. UNVERIFIED / needs a local check

- All throughput speedup figures (no authoritative MB/s table was found; benchmark locally).
- 1.17 release contents; exact fix versions of the threaded-EOF hang relative to 1.16.
- CRAM_OPT_DECODE_MD default behaviour per version.
- Whether the current ExternalProject htslib build links libdeflate.
- CRAM from pipe in all cases.
