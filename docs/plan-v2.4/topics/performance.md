# SAMStat 2.4 - Performance (final topic design)

Status: final design for the implementation workflow. Base: the *minimal* proposal. Grafted from the *researcher* proposal: memory cap as a hard prerequisite, the reader-only-reads pipeline split, the stdin pipe-check use case, and the requirement that the README quotes only benchmark-script numbers.

Sources: `research/code-io.md`, `research/perf-htslib.md`, `research/competitor-speed.md`, `research/code-collect.md`. All throughput figures are from code-io's synthetic-data prototype on one 16-core Apple-silicon Mac with a warm cache. Treat them as ratios, not real-data throughput.

## 0. Principle

- Speed comes from a few changes that remove work or code. Performance adds **no plots and at most one option** (`--reference`, see DECISION NEEDED 2).
- The report must be **byte-identical for any `-t`**. The numbers a researcher reads never depend on how the tool was run.
- **The tool never dies on a long read.** Today a single 1 Mb read costs 3.8 GB RSS and a 96 MB JSON. That is a crash bug, not a speed issue, so it is fixed first.

## 1. Verified code facts this design relies on

Checked against `planning/v2.4` (2.3.0):

| Fact | Where |
|---|---|
| Serial `sam_read1` loop. No `hts_set_threads`. Per-base `seq_nt16_str` copy into a `tl_seq`. | `src/htsinterface/htsglue.c:66-160` |
| `sam_open(filename, "r")`: htslib content detection already exists, but `param.c` rejects inputs by extension first (`.sam .bam .fa .fasta .fna .fq .fastq` plus `.gz`). | `htsglue.c:227`, `param.c:255,287-290` |
| Threads are one per **file**: workers pick the next unprocessed file under a mutex. | `samstat.c:45-150` |
| Batch is a fixed `buffer_size = 10000` records (internal, not a CLI option). | `param.c:390` |
| Every plot is resized to `sb->max_len+1` columns, for R1 **and** R2. | `collect.c:92-116` |
| `stat_collection` holds only `uint64` counts, `plot_data` matrices, and config flags (`collect_end`, `norm_len_plot`, `is_aligned`, `is_partial_report`, `config`, `collect_type`). `stat_collection_finalise()` runs once, in `stat_report.c:45`. So a merge by addition before finalise is sound. Every field still needs an audit (see the merge test). | `collect.h:23-75` |
| `-p` cuts a batch to exactly `top_n` primaries. `n_secondary`/`n_supplementary` are counted in the reader for the **whole batch read**, including records after the cutoff. So under `-p` these two counts depend on the batch size. | `samstat.c:180-215`, `htsglue.c:82-87` |
| FASTA/FASTQ: `prepare_fasta_fastq_input` gz-reads the whole file once to sniff, may write a temp copy, then tld reads it again. | `samstat.c:257ff` |
| The bundled htslib is configured `--disable-bz2 --disable-lzma --disable-libcurl --disable-s3 --without-libdeflate`. `SAMSTAT_USE_SYSTEM_HTSLIB` accepts any system `htslib>=1.16`. | `thirdparty/htslib.cmake` |

The last row has two consequences the proposals missed:
- **CRAM:** CRAMs that use bzip2 or lzma codecs cannot be decoded by the bundled build. Real CRAMs from some pipelines use them (UNVERIFIED how common). Enabling CRAM means enabling bz2 and lzma, which adds two build dependencies.
- **System htslib:** a bioconda build would likely link the system htslib. Any feature needing more than 1.16 must raise the `pkg_check_modules` minimum to match.

## 2. Decisions

### P0. Bounded memory (prerequisite, owned by collect and plots)

- Position axes are bounded **while parsing**, not at render. Long reads use the relative-position axis. The binning rule belongs to the collect and plots topics.
- R2 matrices are allocated lazily, only when a READ2 flag is seen.
- The length histogram is log-binned with bases per bin.
- Performance requires this result: **per-thread stats under about 2 MB, independent of read length.** Without it, W workers multiply 1-4 GB arrays by W. P4 must not land before P0.

### P1. Bump the bundled htslib 1.16 to 1.24 (first, in its own commit)

- Gains, per https://github.com/samtools/htslib/blob/develop/NEWS (taken from `perf-htslib.md`, not re-fetched):
  - faster `bgzf_read` (1.21)
  - CRAM decode speed-ups, and CRAM CVE fixes (1.20, 1.23.1)
  - pipe and threaded-EOF fixes, and no default EBI reference fetch (1.22)
  - more cram2bam work on worker threads (1.24)
- Constraint: since 1.23.1 a thread pool may not be set twice. Create one `hts_tpool` and attach it once per file with `hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p)`.
- Raise the `SAMSTAT_USE_SYSTEM_HTSLIB` minimum to the lowest version the code needs. 1.18 if the modbase block lands (MN tag).
- Refresh the tarball, the MD5 and the release-tarball contents, then run the full ctest suite before anything else builds on it.
- Benchmarks in P6 are taken on this library, not on 1.16.
- libdeflate: see DECISION NEEDED 3.

### P2. htslib decompression threads (about 10 lines)

- One shared pool, attached right after `sam_open` (P1 constraint).
- Measured gain is 1.5-1.8x on BAM and SAM (code-io).
- BGZF only. Plain-gzip FASTQ inflate stays serial in htslib (UNVERIFIED for 1.24; confirm in P6). The FASTQ.gz target is set with that in mind.
- Ships alone, with the `-t` help-text fix ("total threads, also used within one file").

### P3. One reader for everything, content detection, stdin

- **FASTA and FASTQ go through htslib** `sam_open`/`sam_read1`. Delete `prepare_fasta_fastq_input` and the tld reader path for reads.
  - This removes the double inflate, the `/tmp` copy, and the unterminated-last-record workaround.
  - Measured by code-io on 1 M reads: 0.57 s to 0.22 s raw, 1.78 s to 0.78 s gzipped.
  - **Protein FASTA is dropped** (DECISION NEEDED 1).
- **Type by content, not extension.**
  - Delete the extension table in `param.c`.
  - Open with `hts_open` and branch on `hts_get_format(fp)->format` (`sam`, `bam`, `cram`, `fastq_format`, `fasta_format`). An unreadable or unknown file fails only that file (exit 1), not the whole run.
  - This also brings `.cram` (P5) and process substitution for free.
- **`-` means stdin.**
  - Outputs are `stdin.samstat.html`, plus the matching `.json` and `_mqc`, written to `-d` or the current directory. This is the only new naming rule.
  - `-` may appear at most once, and only alongside other inputs whose basename differs.
  - With `-p N` on a pipe, SAMStat stops reading after N primaries. The upstream program gets SIGPIPE, as with `head`. Use case: `minimap2 ... | samstat -p 1000000 -` sanity-checks an alignment before hours of cluster time.
- Quality offset stays Phred+33. FASTQ R1/R2 behaviour must match 2.3 goldens. Whether htslib sets READ1/READ2 from `/1` `/2` name suffixes by default is UNVERIFIED, and must be pinned by a test before the tld path is deleted.

### P4. Reader plus W workers, private accumulators, merge by addition

- **Pipeline:**
  - The reader thread does only `sam_read1` into a batch of recycled `bam1_t` slots.
  - Workers do everything else: the `bam1_t`-to-`tl_seq` conversion, `parse_alignment`, `fix_orientation`, and `collect_stats` into a private `stat_collection`.
  - Moving conversion out of the reader is what should let scaling continue past the 4 workers where the prototype stalled (UNVERIFIED).
- **Queue:** one bounded batch queue (free list plus full list, one mutex, two condition variables). No per-field locks, no work stealing.
- **Batches:** at most 4096 records or 8 Mbases, whichever comes first. Recycled slots whose data buffers exceed 1 MB are shrunk. Worst case is about 64 MB per batch in flight, so about 0.3 GB at W=4.
- **Merge:** new `stat_collection_merge(dst, src)` in `collect.c`. Resize `dst` to the larger dimensions, add element-wise, and add every scalar counter. It runs on the main thread after the workers join and before `stat_report` (and so before finalise).
- **`-p` rule (behaviour change):**
  - The reader counts primaries and stops at exactly N.
  - `n_secondary` and `n_supplementary` count only records read before the cutoff.
  - This makes the counts independent of batch size and thread count. It changes those two numbers in `-p` goldens once, in a dedicated golden-regeneration commit.
- **Gate:** build P4 only if, after P1-P3, a profile of a real 5-10 GB BAM shows per-record work above about 50% of wall time. The prototype (`collect_stats` 37-60% of worker time) says it will.
- **Prototype numbers** at W=4, without merge, report or `-p` (add 0.2-0.6 s for those):

  | Input | Serial | W=4 | Speedup |
  |---|---|---|---|
  | Short BAM | 7.3 s | 1.6 s | 4.6x |
  | Short SAM | 5.6 s | 1.46 s | 3.9x |
  | Long BAM | 6.7 s | 1.38 s | 4.9x |

### P5. Thread budget (one rule, no new options)

- `-t T` means total threads. The default stays 4. There is no core auto-detection: affinity and cgroup detection is unreliable on HPC, and a fixed 4 is predictable.
- With `n_files >= T`, each file runs single-threaded on its own thread, as today. Otherwise files run in sequence, each with `T / n_files` workers.
- The htslib pool shares the budget: `min(4, T)` decompression threads. mosdepth sees no gain past about 4 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6030888).
- `-t 1` runs inline with no thread creation. This is today's code path, so ASAN and debugging stay simple.

### P6. CRAM

- Accept CRAM through P3.
- Reference from `--reference FILE` (`hts_set_fai_filename`) or `REF_PATH`/`REF_CACHE` (DECISION NEEDED 2).
- With a reference available, set `CRAM_OPT_DECODE_MD` so the mismatch profile works. Its default per htslib version is UNVERIFIED, so set it explicitly.
- If the reference is missing, exit 1 with one line naming `--reference` and `REF_PATH`. No silent network fetch: 1.22+ removes the EBI default and the build has `--disable-libcurl`.
- Requires dropping `--disable-bz2 --disable-lzma` from the bundled build (DECISION NEEDED 4).
- What the report shows for MD-less BAMs is owned by the parser and long-read topics, not here.

### P7. Benchmark script

- `scripts/benchmark.sh` (about 30 lines) generates seeded synthetic data with `uv run python`, runs `-t 1/4/8`, and prints reads/s, Mbase/s, wall time and peak RSS.
- The README and the paper quote **only** numbers this script, or a run on a named public dataset, reproduces. This answers the unbenchmarked "2x faster" claim in the paper.

## 3. Targets

All targets are at `-t 4` on the reference Mac, proposals only, and must be re-measured on a real Illumina BAM and a real Dorado BAM before they are frozen or quoted.

| Input | 2.3.0 (measured) | 2.4 target |
|---|---|---|
| Short-read BAM | 0.7 M reads/s | at least 2.5 M reads/s |
| Long-read BAM | 93 Mbase/s | at least 200 Mbase/s (cramino claims about 195 Mbase/s, vendor figure, https://biostars.org/p/9544047, UNVERIFIED) |
| 30x WGS BAM (about 100 Gbase) | about 16-18 min (extrapolated) | 15 min or less (mosdepth parity, PMC6030888) |
| FASTQ, raw | 1.8 M reads/s | at least 2 M reads/s |
| FASTQ.gz (plain gzip) | 0.56 M reads/s | at least 1 M reads/s (serial inflate bound; Falco is about 3x FastQC on gz, https://pmc.ncbi.nlm.nih.gov/articles/PMC7845152) |
| Peak RSS | 3.8 GB on one 1 Mb read | under 1 GB at `-t 4`, independent of read length |
| Determinism | - | byte-identical output for `-t 1/2/4/8` |

## 4. Tests

**Always on (default ctest):**

1. **Merge unit test.** Two halves of a fixture merged must equal the whole, field by field. This is the guard against a non-additive field.
2. **Thread equality.** Extend the existing thread-count test (commit 4824d41): HTML, JSON and `_mqc` byte-identical for `-t 1/2/8` on every golden, plus a generated BAM of more than 100 k reads, so that many batches are in flight.
3. **`-p` boundaries:**
   - N not a multiple of 4096
   - N larger than the file
   - secondary/supplementary counts up to the cutoff only
4. **Stdin:**
   - BAM, SAM, plain FASTQ and gzipped FASTQ piped in equal the file-path result, apart from the output name.
   - Empty stdin gives exit 0 and an empty report.
   - A truncated BAM gives exit 1.
5. **FASTA/FASTQ via htslib.** The existing `fastx` and `golden` labels pass unchanged before the tld path is deleted. Add cases for:
   - an unterminated last record
   - multi-line FASTA
   - `/1` `/2` names
   - a mis-named file (`.txt` holding FASTQ) is accepted
6. **CRAM:**
   - with a reference, counts equal the same data as BAM
   - without one, exit 1 and the named message
   - one CRAM written with lzma or bzip2 codecs
7. **Memory and size bound.** A 1 Mb-read fixture in `tests/longread/` (not `dev/`, which is globbed into four test families) must give RSS under 1 GB and JSON under a fixed ceiling. RSS comes from `/usr/bin/time -l` on macOS and `-v` on Linux. This test may live in the long-read topic; whichever lands first owns it.

**CI:**
- **TSAN job** on ubuntu with `-t 4`. ASAN does not find races.
- The existing ASAN job also runs `-t 4`.

**Perf label** (`-DSAMSTAT_PERF_TESTS=ON`, nightly or manual job only). Gate on ratios, never absolute times:
- `-t 4` time is at most 0.5 of `-t 1` time, on about 2 M generated short reads
- 4x the data takes at most 4.6x the time (linear scaling)
- the long-read fixture with one 500 kb read stays under the RSS ceiling and takes under 3x the time without that read

## 5. Explicitly out

- Region or index parallelism (`.bai`/`.crai`): it misses unmapped reads.
- mmap, SIMD, GPU, a custom BGZF reader.
- A configurable batch size, or a worker count separate from `-t`.
- Core auto-detection, `--fast`, or `--sample-rate`. `-p` covers sampling.
- Progress bars and ETA. One final summary line (reads, seconds, threads) is owned by the CLI topic.
- Duplication estimation.
- Documented remote (S3/HTTP) input.
- A second pass over the data, and in-memory caches.
- Parsing directly from `bam1_t`, to skip the `tl_seq` copy. This is deferred to after 2.4, and only if a profile after P4 shows conversion as the bottleneck.
- Windows.

## 6. Affected files

| File | Change |
|---|---|
| `thirdparty/htslib-1.24.tar.bz2`, `thirdparty/htslib.cmake`, `CMakeLists.txt` | Bump to 1.24. Configure flags (bz2, lzma, libdeflate). Raise the system-htslib minimum. |
| `src/htsinterface/htsglue.c` | Content-sniffing open, shared pool, reader that fills only `bam1_t` batches, conversion function called by workers, CRAM reference and `DECODE_MD`, FASTQ via `sam_read1`. |
| `src/samstat.c` | Thread budget, reader/worker pipeline, merge, `-p` at the reader. Delete `prepare_fasta_fastq_input` and `process_fasta_fastq_file`. |
| `src/collect/collect.c`, `collect.h` | `stat_collection_merge()`. Bounded axes come from P0. |
| `src/param/param.c` | Remove extension detection. Accept `-`. `--reference` (if decided). `-t` help text. `buffer_size` becomes the 4096/8 Mbase rule. |
| `src/report/*` | `stdin` output basename. |
| `tests/cli/`, `tests/fastx/`, `tests/longread/`, `tests/perf/`, `.github/workflows/`, `scripts/benchmark.sh` | Tests and CI as in section 4. |

## 7. Order and effort

Ship each step as its own commit with the suite green.

| # | Step | Effort |
|---|---|---|
| 1 | P1 htslib 1.24 bump, plus regression triage | 1 d |
| 2 | P2 decompression threads, `-t` help | 0.5 d |
| 3 | P3 htslib FASTA/FASTQ, content detection, stdin | 2 d |
| 4 | P0 bounded axes (other topics; hard dependency for step 5) | - |
| 5 | P4 + P5 pipeline, merge, budget, `-p` rule, golden regeneration commit | 3 d |
| 6 | P6 CRAM | 1 d |
| 7 | P7 benchmark script, perf label, TSAN job | 1.5 d |

The total is about 9 days, excluding P0. Steps 1-3 give a 1.5-1.8x gain on BAM, about 2x on FASTQ.gz, plus stdin and CRAM. They also leave less code than 2.3.0.

## 8. Risks

- **Non-additive field:** breaks determinism. Guarded by the merge unit test and the equality test.
- **htslib 1.16 to 1.24 regressions:** isolated as step 1, against 314 tests.
- **htslib FASTQ differences** (R1/R2 flags, comments, quality offset): pinned by tests before the tld path is deleted.
- **Races:** TSAN job.
- **Memory with W workers:** bounded only if P0 lands. Batches are bounded by bases.
- **Targets unverified:** all figures are synthetic, from one machine.

## 9. DECISION NEEDED

1. **Drop protein FASTA?** htslib stores nucleotides only. Recommendation: drop it and note it in the changelog. If you keep it, the tld reader and its two workarounds stay for FASTA only.
2. **How the CRAM reference is supplied.** Options:
   - a `--reference FILE` option, the only new option perf adds
   - `REF_PATH`/`REF_CACHE` environment variables only, with zero new options

   Recommendation: `--reference`, because CRAM is unusable without a reference and environment variables are less discoverable.
3. **libdeflate.** Remove `--without-libdeflate` so htslib uses it when it is found (faster BGZF inflate, size UNVERIFIED)? The cost is builds that differ by machine and one more bioconda dependency. Recommendation: enable it in the bioconda recipe and auto-detect in the bundled build. Output is unaffected.
4. **Enable bz2 and lzma in the bundled htslib** (needed for some CRAMs). This adds libbz2 and liblzma build dependencies. If you say no, CRAM support is "BAM-like CRAMs only", and a codec failure gets a clear error.
5. **`-p` secondary/supplementary semantics change** (counted up to the cutoff only). This regenerates `-p` goldens once. Recommendation: yes.
