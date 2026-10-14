# SAMStat 2.4: work packages for the implementation workflow

Input to a later multi-agent Claude Code run. Specs: `standard-plots.md` (report, data model, JSON v2, MultiQC), `decisions.md` (maintainer rulings), `simplicity-review.md` (budget and cuts). Topic files in `topics/` give background; where they disagree with the three files above, the three files win. Changes made after the plan review are logged in `review-log.md`.

## 0. Rules for every package

1. **File ownership.** A package edits only the files in its "Owns" list. No two packages that run in parallel own the same file. Packages marked **"after X"** run only after X has merged; they may own files X owned. If a change outside the list is needed, the agent stops and reports; it does not edit.
2. **Tests first, green always.** Feature tests are written early as `DISABLED TRUE` (CTest reports "Not Run (Disabled)" and the run stays green). The owning package removes `DISABLED` from its own tests; its gate is that they pass.
3. **ctest is green after every merge**, in both configurations:
   - Release: `cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release && cmake --build build-rel -j && ctest --test-dir build-rel -j4 --no-tests=error --output-on-failure`
   - Sanitizers: same with `-B build-asan -DSAMSTAT_SANITIZE=ON`.
4. **Goldens.** A package never edits `tests/golden/expected/`. A golden diff is allowed only in the two labelled regeneration commits made at gates G3 and G4 ("goldens: regenerate for ..."), via `cmake --build build-rel --target update-goldens`. Any other golden diff rejects the merge.
5. **No scope growth.** Nothing outside `standard-plots.md` and the package goal. No new CLI options (budget: 8 documented options). No new plots.
6. **No CHANGELOG edits** inside packages (conflict magnet). Each package puts one proposed changelog line in its PR description; WP5.4 collects them.
7. Python only for fixture generation, via `uv run python`, never at `ctest` time. Large test inputs are written by CMake scripts (`string(REPEAT)`, `file(WRITE)`) into the build directory.
8. Commits end with the attribution lines from the session's system reminder.

## 1. Agents per package

- **Implementer:** the model in the table. Sonnet for mechanical or well-specified work, Opus where judgement or subtle correctness is involved.
- **Reviewer:** one Opus agent only for packages sized M or L, or with an Opus implementer (marked **[R]** below). It reads the diff against the spec and the "Owns" list, runs the package's tests, and returns PASS or a list of defects. It does not fix. For the other (S, Sonnet) packages the gate verifier's ctest run plus the golden rule is the review.
- **Gate verifier:** one Sonnet agent per automated gate (G1-G4) does the merges in the stated order, runs both ctest configurations after each merge, and checks the golden rule. G0 and G5 are done by the maintainer.

23 implementers, 15 reviewers and 4 gate verifiers: 42 agent runs. Within a wave, packages without an "after" mark run in parallel.

## 2. Overview

| Wave | Packages | Gate |
|---|---|---|
| 0 | WP0.1 release 2.3.0 prep · WP0.2 tld strbuf fix upstream | G0: maintainer tags v2.3.0, opens the bioconda PR, pushes the tld fix |
| 1 | WP1.1 test infra · WP1.2 fixtures · WP1.3 htslib bump · WP1.4 CLI cleanup · WP1.5 tld pin + JSON flush removal · WP1.6 SVG primitives · WP1.7 report text | G1: goldens byte-identical |
| 2 | WP2.1 pending tests · WP2.2 parser (identity, X ops, NM reset) · WP2.3 bench script + TSAN job | G2: goldens byte-identical |
| 3 | WP3.1a collector v2 · WP3.2 htslib reader, content detection, stdin; then WP3.1b JSON v2 preview (after 3.1a) | G3: golden regeneration #1 |
| 4 | WP4.1 intra-file threading · WP4.2a report page · WP4.2b plot renderers · WP4.3 JSON v2 + MultiQC | G4: golden regeneration #2 (format_version 2) |
| 5 | WP5.2 CRAM · WP5.3 stderr line + progress · WP5.4 docs + changelog; then WP5.1 dead-code cleanup (after 5.3) | G5: release gate, maintainer tags v2.4.0 |

Critical path: WP1.2 -> WP2.1 -> WP3.1a -> WP3.1b -> WP4.2b/4.2a -> WP5.3 -> WP5.1 -> G5.

**If the schedule slips, WP5.2 (CRAM) is dropped first** and moves to 2.5; nothing else depends on it.

## 3. File ownership by wave (conflict check)

`a/b` in a cell means sequential packages in that wave (b after a).

| File / directory | W0 | W1 | W2 | W3 | W4 | W5 |
|---|---|---|---|---|---|---|
| `CHANGELOG.md`, `ChangeLog`, `NEWS` | 0.1 | | | | | 5.4 |
| `packaging/bioconda/meta.yaml` | 0.1 | | | | | 5.2 |
| `thirdparty/htslib-*.tar.bz2`, `thirdparty/htslib.cmake` | | 1.3 | | | | 5.2 |
| `src/tld` (pin), `tests/cmake/check_tld_pin.cmake` | | 1.5 | | | | |
| `src/report/json_report.c` | | 1.5 | | 3.1b | 4.3 | |
| `src/param/param.{c,h}` | | 1.4 | | 3.2 | 4.1 | 5.1 |
| `Readme.org` | | 1.4 | | | | 5.4 |
| `src/plot/svg_plot.{c,h}`, `src/plot/svg_plot_test.c`, `tests/svg/` | | 1.6 | | | 4.2b | |
| `src/plot/report_plots.{c,h}`, `src/plot/report_plots_test.c` (new) | | | | | 4.2b | |
| `src/plot/CMakeLists.txt` | | 1.6 | | | 4.2b | 5.1 |
| `src/plot/plot.{c,h}`, `plot_group.*`, `plot_interval.*`, `plot_test.c` | | | | | | 5.1 |
| `src/report/report_help.{c,h}`, `src/CMakeLists.txt` | | 1.7 | | | 4.2a (`report_help.*` only) | |
| `tests/cmake/check_json.cmake`, `compare_golden.cmake` | | 1.1 | | | | |
| `tests/golden/CMakeLists.txt` | | 1.1 | 2.1 | 3.1b | | |
| `tests/json/CMakeLists.txt`, `tests/json/golden/` | | 1.1 | | | 4.3 | |
| `tests/packaging/CMakeLists.txt` | | 1.1 | | | | |
| `tests/tools/make_longread_fixtures.py`, `tests/data/longread/` | | 1.2 | | | | |
| `tests/cli/CMakeLists.txt` | | 1.4 | 2.1 | 3.2 | 4.1 | 5.3 |
| `tests/longread/CMakeLists.txt`, `tests/cmake/run_ultralong.cmake` | | | 2.1 | 3.1b | 4.3 | |
| `tests/cmake/run_threads.cmake` | | | 2.1 | | 4.1 | |
| `tests/report/CMakeLists.txt`, `tests/cmake/check_plot_set.cmake`, `tests/cmake/run_html_size.cmake` | | | 2.1 | | 4.2a | |
| `tests/reader/CMakeLists.txt` | | | 2.1 | 3.2 | 4.3 | 5.2 |
| `tests/fastx/CMakeLists.txt` | | | 2.1 | 3.2 | 4.3 | |
| `src/sambamparse/sam_bam_parse.{c,h}`, `tests/reader/unit_parse_alignment.c` | | | 2.2 | | | |
| `src/htsinterface/htsglue.{c,h}` | | | 2.2 | 3.2 | 4.1 | 5.2 |
| `tests/perf/bench.sh`, `docs/benchmark.md` | | | 2.3 | | | |
| `.github/workflows/cmake.yml` | | | 2.3 | | | 5.2 |
| `src/collect/collect.{c,h}`, `src/convert_tables.c` | | | | 3.1a | | 5.1 |
| `tests/collect/CMakeLists.txt`, `tests/collect/unit_collect.c` (new) | | | | 3.1a/3.1b | 4.3 | |
| `docs/samstat.schema.json`, `docs/json-format.md` | | | | 3.1b | 4.3 | |
| `tests/cmake/run_golden.cmake`, `run_expect.cmake` | | | | 3.1b | | |
| `tests/cmake/run_json.cmake`, `run_collect.cmake`, `check_schema_keys.cmake` | | | | 3.1b | 4.3 | |
| `src/samstat.c` | | | | 3.2 | 4.1 | 5.3/5.1 |
| `src/thread/thread_data.{c,h}` | | | | | 4.1 | |
| `src/report/stat_report.{c,h}` | | | | 3.2 (one hunk) | 4.2a | |
| `src/report/mqc_report.c`, `docs/multiqc.md`, `tests/multiqc/CMakeLists.txt`, `tests/cmake/run_mqc.cmake` | | | | | 4.3 | |
| `tests/data/cram/` | | | | | | 5.2 |
| `tests/golden/expected/` | | | | gate G3 | gate G4 | |

---

## 4. Packages

### Wave 0 (before any code change)

#### WP0.1 Release 2.3.0 prep. Sonnet, S
- **Goal:** make 2.3.0 releasable as it is (release-cleanup R1, R2, R8).
- **Owns:** `CHANGELOG.md`, `ChangeLog`, `NEWS`, `packaging/bioconda/meta.yaml`.
- **Do:** date `## [2.3.0]`; add "Superseded by CHANGELOG.md" at the top of `ChangeLog` and `NEWS`; in `meta.yaml` drop host deps `bzip2`, `xz`, `libcurl` (the bundled htslib is built without them, `thirdparty/htslib.cmake`), remove the DRAFT header, add a test command that runs `samstat` on a 3-read FASTQ written with `printf`. The sha256 stays a placeholder until the release asset exists.
- **Acceptance:** `packaging:metadata` passes; offline tarball check: `packaging/make_source_tarball.sh 2.3.0 dist`, extract into an empty temp dir, build with networking off, `samstat -v` prints `SAMStat 2.3.0`, `ctest -L "build|cli"` green.
- **Tests:** none new.

#### WP0.2 tld strbuf fix (upstream repo). Sonnet, S
- **Goal:** fix `tld_strbuf_resize` never updating `alloc_len` (release-cleanup R3).
- **Owns:** a patch against https://github.com/timolassmann/tld only (prepared in a scratch clone, not in `src/tld`).
- **Do:** add `b->alloc_len = new_size;`; add an upstream test doing 1e6 appends that asserts `alloc_len >= len` and fewer than 64 reallocations.
- **Acceptance:** upstream test passes. The maintainer reviews and pushes (agents never push).

#### Gate G0 (maintainer)
- Sign off `decisions.md`.
- Tag `v2.3.0` on main; the release workflow attaches tarball and `.sha256`. Put the real sha256 in `meta.yaml` and open the bioconda PR (https://github.com/bioconda/bioconda-recipes).
- Push the WP0.2 tld commit upstream; note its hash for WP1.5.
- Create the integration branch `develop/v2.4` from the v2.3.0 tag plus this planning directory.

### Wave 1 (foundation; seven packages in parallel)

#### WP1.1 Test infrastructure. Sonnet, M [R]
- **Goal:** testing T1.
- **Owns:** `tests/cmake/check_json.cmake`, `tests/cmake/compare_golden.cmake`, `tests/golden/CMakeLists.txt`, `tests/json/CMakeLists.txt`, `tests/json/golden/` (deleted), `tests/packaging/CMakeLists.txt`. Moves the 4 files of `tests/json/golden/` into `tests/golden/expected/` (a move, not a content change; the only exception to rule 4).
- **Do:** `>=` and `<=` operators in `check_json`; `compare_golden` reports the first differing dotted path; one golden mechanism (`update-goldens`); fixture budget check (fails if `tests/data` > 1 MB).
- **Acceptance:** all 314 existing tests still pass; a deliberately broken golden prints its path.
- **Tests:** `packaging:fixture_budget`; `json:check_json_ops` (self-test of `>=`/`<=` on `tests/json/data/mini.json`).

#### WP1.2 Long-read fixtures and generator. Opus, M [R]
- **Goal:** testing T2, ground-truth fixtures.
- **Owns:** `tests/tools/make_longread_fixtures.py`, `tests/data/longread/` (incl. `README.md` with expected values).
- **Do:** stdlib-only, seeded generator; `--reads N` option for `bench.sh` (writes SAM). Files: `ont_like.sam` + `.bam` (about 40 reads, 500 bp-30 kb, about 8% error, mixed MAPQ, soft clips, MD and NM; the BAM adds one read with > 65535 CIGAR ops), `supp_chain.sam`, `hifi_like.sam` (Q up to 60), `long.fq` (unterminated last record; two reads with equal Phred mean but different error-probability mean), `eqx_nomd.sam` (= and X CIGAR ops, no MD), `nonm.sam` (no NM), `ubam_like.sam` (no `@SQ`, every record flag 4), `sq_nomap.sam` (`@SQ` present, every record flag 4), `mapq255.sam` (mapped reads with MAPQ 255), `single_len.sam` (all reads 150 bp, paired). Plus hand-computable cases (N50 of 5 reads, identity with one indel, soft-clip percentage). BAMs made with htslib `test/test_view` as in `tests/data/README.md`.
- **Acceptance:** total <= 300 KB, each file <= 100 KB; expected values come from generator ground truth, never from samstat; `samtools view` (if available) parses every file.
- **Tests:** none registered (WP2.1 uses them).

#### WP1.3 htslib bump. Sonnet, S
- **Goal:** bundle the latest htslib release (target 1.24; release date per `research/perf-htslib.md`, UNVERIFIED; check https://github.com/samtools/htslib/releases), behaviour unchanged.
- **Owns:** `thirdparty/htslib-<ver>.tar.bz2` (replaces 1.16), `thirdparty/htslib.cmake`.
- **Do:** switch `URL_MD5` to `URL_HASH SHA256=`; keep configure flags unchanged; system-htslib minimum stays `>=1.16`.
- **Acceptance:** full suite green, goldens byte-identical. Any golden diff is explained in the PR and goes to the maintainer before merge.

#### WP1.4 CLI cleanup. Sonnet, M [R]
- **Goal:** cli-ux decisions with simplicity-review C13, C14, C16.
- **Owns:** `src/param/param.{c,h}`, `Readme.org`, `tests/cli/CMakeLists.txt`.
- **Do:**
  - `-l`, `--seed`, `--norm`, `--plotend` accepted and ignored, each printing `samstat: warning: <opt> is deprecated and ignored`, removed from `-h`. Fixes the 2.3 `--seed` abort (no `case` for `OPT_SEED`, rc 134). `--plotend`/`--norm` set nothing (so no `_end`/`len_norm` matrices are allocated); the `collect_end`/`norm_len_plot` fields stay until WP5.1.
  - `--threads` long name, `--nthreads` kept as an alias. Help text unchanged in meaning until WP4.1.
  - `-d` creates the directory (`mkdir -p`) instead of erroring (`param.c:199-208`).
  - Up-front output checks before any record is read: unwritable output dir fails that input (exit 1) with `samstat: error: cannot write <path>: <strerror>; use -d DIR to write elsewhere.`; two inputs resolving to the same output path is exit 2, naming both.
  - New `print_help()` (the 8 options of `standard-plots.md`/cli-ux §7, `-d` instead of `-o`); `Readme.org` matches.
- **Acceptance:** exit codes 0/1/2 unchanged in meaning; goldens byte-identical.
- **Tests (label `cli`):** `deprecated_l`, `deprecated_seed`, `deprecated_norm`, `deprecated_plotend` (exit 0, one warning line, JSON identical to a run without the flag); `dir_created`; `outdir_unwritable` (exit 1, within 1 s, message contains `-d DIR`); `input_dir_readonly` (input copied into a `chmod 555` build-dir directory, no `-d`: exit 1, message contains `-d DIR`; skipped when run as root); `output_collision` (exit 2); `version_regex` (`^SAMStat [0-9]+\.[0-9]+\.[0-9]+$`); `help_matches_readme` (every long option in `-h` appears in `Readme.org`).

#### WP1.5 tld pin and JSON flush removal. Sonnet, S
- **Depends:** G0 (tld fix pushed).
- **Owns:** `src/tld` (submodule pin), `tests/cmake/check_tld_pin.cmake`, `src/report/json_report.c`.
- **Do:** bump the pin once; delete `JSON_BUF_BYTES` and the periodic `json_flush` calls (`json_report.c:31`); build the whole JSON in one strbuf, write once.
- **Acceptance:** JSON byte-identical for every golden.
- **Tests:** `packaging:tld_pin` updated.

#### WP1.6 SVG primitives. Opus, M [R]
- **Goal:** presentation D1 primitives, standalone.
- **Owns:** `src/plot/svg_plot.{c,h}` (new), `src/plot/svg_plot_test.c` (new), `src/plot/CMakeLists.txt`, `tests/svg/CMakeLists.txt` (new area, auto-registered by `tests/CMakeLists.txt:36`).
- **Do:** about 350 lines: nice-ticks; linear and log-10 axes (1-2-5 ticks); one `<path>` per series with dash pattern and gaps for missing values, one series-level `<title>`; step histogram; vertical marker line with label; inline swatch row (for plot 4 and the MAPQ table); panel grid (rows x columns with shared x, optional y-scale shared per row); `role="img"`, `<title>`, `<desc>`. Output into a `tld_strbuf`. No JavaScript, no per-point elements, no direct line labels.
- **Acceptance:** unit test renders each primitive; deterministic bytes for fixed input; output is well formed (checked by a small tag-balance and attribute-quoting checker inside `svg_plot_test.c`; no external XML tool).
- **Tests (label `report`):** `svg:unit` (runs `svg_plot_test`, including the well-formedness check).

#### WP1.7 Report text. Opus, S [R]
- **Goal:** single source of all report wording.
- **Owns:** `src/report/report_help.{c,h}` (new), `src/CMakeLists.txt` (add `report/report_help.c` to `source_files`).
- **Do:** a static table indexed by plot id: title, caption, three "How to read this" lines, N/A texts, note templates, MAPQ labels and the MAPQ table caption; text exactly as `standard-plots.md` sections 1.3, 1.4, 4 and 5.
- **Acceptance:** compiles and links; maintainer reads the text once (it is user-facing).
- **Tests:** none (WP4.2a tests the rendered text).

#### Gate G1
Merge order: 1.1, 1.3, 1.5, 1.4, 1.2, 1.6, 1.7. After each: both ctest configurations green, goldens byte-identical (1.1's move excepted).

### Wave 2 (tests ahead of features; parser)

#### WP2.1 Pending feature tests. Sonnet, M [R]
- **Depends:** WP1.1, WP1.2.
- **Owns:** `tests/longread/CMakeLists.txt` (new), `tests/cmake/run_ultralong.cmake` (new), `tests/cmake/run_threads.cmake`, `tests/cli/CMakeLists.txt`, `tests/report/CMakeLists.txt`, `tests/cmake/check_plot_set.cmake` (new), `tests/cmake/run_html_size.cmake` (new), `tests/reader/CMakeLists.txt`, `tests/fastx/CMakeLists.txt`, `tests/golden/CMakeLists.txt`.
- **Do (all feature tests `DISABLED TRUE`, with a comment naming the package that enables them):**
  - `longread:*` value tests against `standard-plots.md` section 7 keys under the transitional `v2.` prefix (testing T3): exact counts; N50 and median within 1%; median identity within 0.5 pp; mean Q within 0.1 on `long.fq`; MAPQ bins sum to primary reads; `pct_clipped` on the hand-computable case; `ubam_like.sam` gives `is_aligned 0`, `sq_nomap.sam` gives `is_aligned 1`; `mapq255.sam` counts under `mapq_ge30`. -> WP3.1b
  - `longread:ultralong` (T4): `run_ultralong.cmake` writes a 1 Mb single-read SAM into the build dir with `string(REPEAT)`; asserts JSON < 1 MB, every absolute-axis array length <= 1000, every relative-axis array <= 100, identity arrays <= 1001. -> WP3.1b
  - `report:ultralong_html` (HTML < 1 MB on the 1 Mb read) and `report:short_html_size` (`run_html_size.cmake` writes a paired aligned SAM of 999 bp reads, enough records to fill every bin, HTML < 1 MB). -> WP4.2a
  - `cli:threads_one_file` (T6): `run_threads.cmake` runs `-t 1/3/8` with `SAMSTAT_BATCH_RECORDS=8` and `SAMSTAT_BATCH_BASES=20000` (test-only, section 4 WP3.2) on `paired_small.sam` and `ont_like.bam` (at least 5 batches each), byte-identical JSON; asserts the `--verbose` line `workers: N` shows N > 1 for `-t 3` and `-t 8`. -> WP4.1
  - `report:plot_set` (T5): titles in order for a short BAM, `ont_like.bam`, `long.fq` (plots 1-4 only); no `NaN|Infinity|undefined`; short-layout axis ends at 150 for `single_len.sam`; `report:offline`: no `https?://` in `src=`/`href=`. -> WP4.2a
  - `cli:stdin_*`, `cli:unknown_input_per_file`, `fastx:sniff_misnamed`. -> WP3.2
  - Enabled now (pass on 2.3): `reader:truncation_loop` (cut `ont_like.bam` at 5 offsets; exit 1 or 2, no signal, no report left); `fastx:r1r2_suffix` (FASTQ with names ending `/1` and `/2`: `read2_reads` = 0, pins 2.3 behaviour); missing `fastx` edge cases (no trailing newline plain and `.gz`, empty `.gz`, multi-line FASTA) where not already covered.
  - Golden glob extended to `tests/data/longread/*.{sam,bam,fq}`; their goldens are created at G2 by the gate verifier (new files, not a diff).
- **Acceptance:** `ctest` green; disabled count printed in the PR.

#### WP2.2 Parser: identity fields, X ops, NM reset. Opus, M [R]
- **Owns:** `src/sambamparse/sam_bam_parse.{c,h}`, `src/htsinterface/htsglue.{c,h}`, `tests/reader/unit_parse_alignment.c`.
- **Do:** (`standard-plots.md` 3.4, 3.5)
  - In the first CIGAR loop of `parse_alignment`: per-record `M`, `I_ev`, `I_bases`, `D_ev`, `D_bases`, `S_bases`, `X_ops` into `struct aln_data` (`htsglue.h`).
  - X ops mark mismatches when there is no MD; per-record flag `has_md_or_x`.
  - `read_bam_chunk`: reset `a->error = -1` before the NM lookup.
  - The collector does not consume the new fields yet.
- **Acceptance:** goldens byte-identical (no fixture uses =/X ops today; checked).
- **Tests (label `reader`):** cases in `unit_parse_alignment.c`: identity counts for a CIGAR with one insertion and one deletion; soft-clip count; `=`/X CIGAR without MD gives mismatches; a record without NM after one with NM has `error == -1`.

#### WP2.3 Benchmark script and TSAN job. Sonnet, S
- **Depends:** WP1.2.
- **Owns:** `tests/perf/bench.sh` (new), `docs/benchmark.md` (new), `.github/workflows/cmake.yml`.
- **Do:**
  - POSIX sh script (testing T7). Generates **SAM** input in the build dir via the generator (about 2 M short reads, about 1 Gb ONT-like reads, and 8 reads of 1 Mb); records generation time separately from the timed runs. Runs `-t 1/4/8` under `/usr/bin/time` (macOS `-l`, Linux `-v`); prints wall time, peak RSS, Mreads/s, Mbases/s, the htslib version and whether it was built with libdeflate. Also times 2.3.0 (a path given as an argument) on the same files.
  - `docs/benchmark.md`: a results table for the maintainer, with columns for machine, OS, build (bundled or system htslib, libdeflate yes/no). Targets are ratios, measured on x86-64 Linux after WP3.1a: `-t 4` >= 2.0x and `-t 8` >= 3.5x the `-t 1` rate of the same build; FASTQ.gz >= 1.5x the 2.3.0 rate; 8 x 1 Mb reads at `-t 8` < 1 GB peak RSS. The prototype numbers in `topics/performance.md` are not quoted.
  - Add an `ubuntu-tsan` CI job (`-fsanitize=thread`) that runs only `ctest -L golden` and `ctest -R threads`, with a step `sudo sysctl vm.mmap_rnd_bits=28` before the tests (TSAN aborts with "unexpected memory mapping" on recent runner kernels otherwise, https://github.com/google/sanitizers/issues/1716, UNVERIFIED, not fetched).
- **Acceptance:** script runs on macOS; TSAN job green on 2.3 code. No perf CTest label, no gates.

#### Gate G2
Merge 2.2, 2.1, 2.3. The verifier then creates goldens for the new `tests/data/longread/` inputs (additions only, one commit "goldens: add long-read fixtures"). Both configurations green.

### Wave 3 (data model and reader)

#### WP3.1a Collector v2. Opus, L [R]
- **Depends:** WP2.1, WP2.2, WP1.4.
- **Owns:** `src/collect/collect.{c,h}`, `src/convert_tables.c`, `tests/collect/CMakeLists.txt`, `tests/collect/unit_collect.c` (new).
- **Do:** (`standard-plots.md` section 3)
  - All accumulators of 3.1-3.6 as fixed-size integer arrays, array-of-structs per position bin; `E[q]` table generated by `convert_tables.c`; per-read Q bin by binary search; relative bin by an incremental walk (no division per base).
  - Input contract of 3.0: decisions per record; delete the `config_stats` first-chunk `collect_type` logic (`collect.c:248-276`) now, not in wave 5. `is_aligned` is computed in finalise from `has_sq` and mapped counts (3.6); the reader supplies `has_sq`.
  - Legacy 2.3 position matrices capped at 1000 columns and no longer resized to `max_len` (`collect.c:92-116`); legacy len and Q-matrix kept only so the 2.3 report still renders until WP4.2a. `struct plot_data` and `plot_data_alloc/resize_len/free` are used unchanged.
  - `stat_collection_merge(dst, src)` with the per-field semantics of 3.6 (sum, OR, min/max; `min_length` starts at `UINT64_MAX`), including the legacy matrices.
  - Finalise function computing every summary value (3.7).
  - `stat_collection_config_additional_plot` stays (called from `samstat.c`); WP5.1 removes it.
  - Report in the PR: single-thread Mbases/s from `bench.sh -t 1` on the long-read SAM before and after; the v2 collector may cost at most 25% more wall time than 2.3.
- **Acceptance:** `collect:unit` passes; existing tests green; any v1 golden change is limited to legacy position arrays truncated at 1000 and to files whose first chunk lacked qualities (listed in the PR).
- **Tests (label `collect`):** `collect:unit` in `unit_collect.c`: merge of a filled collection with a zero-read collection equals the filled one; merge is order-independent; every storage length edge >= 1000 is also a display edge; `E[0] = 2^30`; identity bin of a hand-computed CIGAR; mean Q from a synthetic collection whose `qual_error_sum_all` exceeds `UINT64_MAX`.

#### WP3.1b JSON v2 preview. Sonnet, M [R] (after 3.1a)
- **Owns:** `src/report/json_report.c`, `docs/samstat.schema.json`, `docs/json-format.md`, `tests/cmake/run_golden.cmake`, `run_json.cmake`, `run_collect.cmake`, `run_expect.cmake`, `check_schema_keys.cmake`, `tests/collect/CMakeLists.txt`, `tests/golden/CMakeLists.txt`, `tests/longread/CMakeLists.txt`, `tests/cmake/run_ultralong.cmake`.
- **Do:** emit the complete v2 layout (section 7) under a top-level `"v2"` key; v1 keys unchanged. Schema and `json-format.md` document `"v2"` as a transitional optional property. Remove the 4 MB golden skips (`tests/golden/CMakeLists.txt:12`, `run_golden.cmake`). Enable the `longread:*` tests, including `longread:ultralong`.
- **Acceptance:** enabled tests pass; expected golden diff = the added `"v2"` object plus the 3.1a diffs. Nothing else.
- **Tests (label `collect`, `longread`):** as enabled; plus `collect:mapq_sum` (MAPQ class counts sum to primary reads) and `collect:axis_rule` (150 bp median gives `layout short`, 10 kb gives `long`; both axes present).

#### WP3.2 One htslib reader, content detection, stdin. Opus, L [R]
- **Depends:** WP1.3, WP1.4, WP2.1.
- **Owns:** `src/htsinterface/htsglue.{c,h}`, `src/samstat.c`, `src/param/param.{c,h}`, `tests/cli/CMakeLists.txt`, `tests/fastx/CMakeLists.txt`, `tests/reader/CMakeLists.txt`, and one hunk of `src/report/stat_report.c` (skip `samstat_file_stat_section`, `stat_report.c:52`, when the input is `-`; WP4.2a rewrites the file later).
- **Do:** (performance P2, P3, `-p` rule from P4; `standard-plots.md` 3.0)
  - FASTA/FASTQ through `sam_open`/`sam_read1`; delete `prepare_fasta_fastq_input` (`samstat.c:257`) and the tld read path. When `hts_get_format(fp)->format` is `fasta_format` or `fastq_format`: no `aln_data`, `parse_alignment` not called, `BAM_FPAIRED|BAM_FREAD1|BAM_FREAD2` ignored, every read is R1 (htslib sets these flags from `/1`, `/2` name suffixes in `fastq_parse1`, htslib 1.16 `sam.c`). Remove `s->is_aligned = 1` (`samstat.c:178`); set `has_sq` from the header.
  - Empty input rule (3.0): zero-byte or whitespace-only input after decompression is exit 0 with a warning and `total_reads 0`, checked before format detection; empty-sequence records are not counted. Protein FASTA is read as nucleotides (non-ACGT letters become N); the changelog line says it is unsupported.
  - Type from content (`hts_get_format`), extension table in `param.c` deleted; an unknown non-empty input fails only that input (exit 1). Update `cli:ext_unknown_*` (rc 2 to rc 1, other report present).
  - `-` = stdin, at most once (second `-` is exit 2). Output stem `stdin` via a new helper in `param.c` (`samstat_output_stem(p, id)`); the report writers switch to it in WP4.2a/WP4.3.
  - One shared `hts_tpool` attached once per file (htslib >= 1.23.1 forbids setting it twice, per topic research, UNVERIFIED), size per the WP4.1 formula (1 thread at the default until WP4.1 lands).
  - Batch = 4096 records or 8 Mbases, whichever first (`param.c:390` today 10000 records). The two limits are read in one place and can be lowered by the undocumented test-only environment variables `SAMSTAT_BATCH_RECORDS` and `SAMSTAT_BATCH_BASES` (not a user option; not in `-h` or docs).
  - `-p N`: the reader stops mid-batch at exactly N primaries; `secondary`/`supplementary` are counted only in the reader and only for records read before the cutoff.
- **Acceptance:** `fastx`, `golden` labels unchanged except `-p` goldens' secondary/supplementary counts (if any exist); `fastx:json:newline_only`, `header_only`, `empty_gz` and `cli:empty_fq` still pass unchanged; `fastx:r1r2_suffix` still passes; stdin tests enabled (with `-d`, since report names switch in wave 4: assert on JSON content only).
- **Tests:** enable `cli:stdin_*`, `cli:unknown_input_per_file`, `fastx:sniff_misnamed`; add `cli:stdin_double_dash` (exit 2), `reader:p_cutoff_counts` (N not a multiple of 4096; N > file size).

#### Gate G3 (golden regeneration #1)
Merge 3.2 (expect no golden diff, or only `-p` secondary/supplementary counts), then 3.1a, then 3.1b. Then one commit "goldens: regenerate for v2 preview and -p cutoff": the verifier runs `update-goldens` and checks the diff is limited to the expectations of the three packages. Both configurations green.

### Wave 4 (threading, report, JSON v2)

#### WP4.1 Intra-file threading. Opus, L [R]
- **Depends:** WP3.1a (merge), WP3.2.
- **Owns:** `src/samstat.c`, `src/thread/thread_data.{c,h}`, `src/htsinterface/htsglue.{c,h}`, `src/param/param.{c,h}`, `tests/cli/CMakeLists.txt`, `tests/cmake/run_threads.cmake`.
- **Do:** (performance P4, P5)
  - Reader thread does `sam_read1` into recycled `bam1_t` batches, counts secondary/supplementary and applies the `-p` cutoff (never the workers). W workers convert, parse, collect into private `stat_collection`s. One bounded queue (one mutex, two condition variables) with a fixed pool of **W+2 batches**; merge after join, before finalise.
  - **Thread budget, `-t T` = total threads**, default 4, no core detection: `T = 1` runs inline, no pool. `T = 2`: reader + 1 worker, no decompression pool. `T >= 3`: decompression pool `P = max(1, T/4)`, workers `W = T - 1 - P` (T=4: P 1, W 2; T=8: P 2, W 5). With `n_files >= T`, files run one per thread, each inline without a pool (2.3 scheduling); otherwise files run one after another, each with the full budget. `--verbose` prints `workers: W, decompression threads: P`. `-t` help text: "total threads (default 4)".
  - Memory bound (write it into `docs/benchmark.md` via the PR description for WP5.4): at most W+2 batches of <= 8 Mbases plus one longest read each, plus W worker buffers sized by the longest read, plus W x 0.85 MB collections. Estimate at `-t 8` for 1 Mb reads: about 7 x (8 + 1) Mb x about 4 B per base = about 250 MB (UNVERIFIED; `bench.sh` measures it).
- **Acceptance:** output byte-identical for any `-t` (JSON, HTML, `_mqc`); goldens unchanged; `reader:p_cutoff_counts` passes at `-t 1` and `-t 8`; TSAN and ASAN jobs green.
- **Tests:** enable `cli:threads_one_file`; extend it to HTML and `_mqc` outputs; add a `-t 8` variant of `reader:p_cutoff_counts`.

#### WP4.2b Plot renderers. Opus, M [R]
- **Depends:** WP3.1a, WP1.6.
- **Owns:** `src/plot/report_plots.{c,h}` (new), `src/plot/report_plots_test.c` (new), `src/plot/svg_plot.{c,h}`, `src/plot/svg_plot_test.c`, `src/plot/CMakeLists.txt`, `tests/svg/CMakeLists.txt`.
- **Do:** the six plots and the MAPQ swatches of `standard-plots.md` sections 2 and 4, from the finalised collection, on top of the `svg_plot` primitives. Fixed interface (WP4.2a codes against it in parallel):
  - `int report_plot_render(tld_strbuf *o, const struct stat_collection *s, int plot_id);` for plot_id 1-6, writing either the figure or the N/A line;
  - `int report_mapq_swatch(tld_strbuf *o, int mapq_class);`
  - the two palettes as `static const char *` arrays in `report_plots.h`.
  Adds files only; the Plotly code in `plot.c`, `plot_group.*`, `plot_interval.*` stays until WP5.1.
- **Acceptance:** deterministic bytes; no per-point elements; x ranges per 3.2; plot 5 rows share y between R1 and R2.
- **Tests (label `report`):** `svg:report_plots` (runs `report_plots_test`: each plot renders for a synthetic collection; the base palette and MAPQ palette share no hex value).

#### WP4.2a Report page. Opus, M [R]
- **Depends:** WP3.1a, WP1.7; codes against the WP4.2b interface.
- **Owns:** `src/report/stat_report.{c,h}`, `src/report/report_help.{c,h}`, `tests/report/CMakeLists.txt`, `tests/cmake/check_plot_set.cmake`, `tests/cmake/run_html_size.cmake`.
- **Do:** the page of `standard-plots.md` section 1 (header, strip, notes, MAPQ table, plot calls, footer without command line, print CSS) from the finalised values. `stat_report.c` stops calling the Plotly emitters, the dropdown, the Plotly MAPQ table (`stat_report.c:547-581`) and the Plotly credit line (`stat_report.c:628`). Report file name via `samstat_output_stem()`.
- **Acceptance:** report for each fixture opens offline; HTML < 1 MB for the 1 Mb read and for the 999 bp paired file; headline numbers consistent with JSON.
- **Tests (label `report`):** enable `report:plot_set`, `report:offline`, `report:ultralong_html`, `report:short_html_size`; add `report:na_notes` (FASTQ and `ubam_like.sam` show plots 1-4 only, no MAPQ table and no plot 5-6 text; `sq_nomap.sam` shows the MAPQ table and `No mapped reads.` for plots 5-6; `eqx_nomd.sam` shows mismatches from X; a no-MD, no-X SAM shows the mismatch N/A line and insertion/deletion rows; `nonm.sam` shows the identity note), `report:svg_a11y` (every svg has `role="img"`, `<title>`, `<desc>`), `report:no_flags` (no `PASS|WARN|FAIL`), `report:partial_header`.

#### WP4.3 JSON v2 and MultiQC. Sonnet, M [R]
- **Depends:** WP3.1b, WP3.2.
- **Owns:** `src/report/json_report.c`, `src/report/mqc_report.c`, `docs/json-format.md`, `docs/samstat.schema.json`, `docs/multiqc.md`, `tests/json/CMakeLists.txt`, `tests/multiqc/CMakeLists.txt`, `tests/cmake/run_mqc.cmake`, `tests/longread/CMakeLists.txt`, `tests/cmake/run_ultralong.cmake`, `tests/collect/CMakeLists.txt`, `tests/reader/CMakeLists.txt`, `tests/fastx/CMakeLists.txt`, `tests/cmake/run_collect.cmake`, `tests/cmake/run_json.cmake`, `tests/cmake/check_schema_keys.cmake`.
- **Do:**
  - Promote `"v2"` to top level, delete v1 keys, `format_version: 2`; summary floats `%.4f`, `null` when not applicable.
  - Rewrite every v1 `CHECKS` path in the owned test files to its v2 equivalent (e.g. `read1.length.*`, `read1.base_quality.*` in `tests/collect`, 27 refs in `tests/reader`, 9 in `tests/fastx`), or delete the check when the quantity no longer exists (Q x position matrix, clip-trimmed length). `run_collect.cmake` and `run_json.cmake` expect `format_version` 2. List every changed or deleted test in the PR.
  - MultiQC per `standard-plots.md` section 8 (7 columns, 5 sections, display-binned length, relative axis for both position linegraphs, `samstat_error`). File names via `samstat_output_stem()`. `fn_clean_sample_names: false` note moved to the top of `docs/multiqc.md`.
  - `json_report.c` and `mqc_report.c` stop reading `struct plot_data`.
- **Acceptance:** schema validates every golden; `check_schema_keys` passes; `longread:*` tests switched from `v2.` paths to top-level paths.
- **Tests:** `json:format_version_2`; `multiqc:general_ids` (exact id list); `multiqc:error_absent_fastq`; `multiqc:relative_axis` (both position sections have 100 x values for a short and a long fixture); hand-computed `samstat_n50`, `samstat_median_identity` on fixtures; `multiqc:run` asserts the new section for all samples (CI multiqc job).

#### Gate G4 (golden regeneration #2)
Merge 4.1 (goldens unchanged), 4.2b, 4.2a, then 4.3, then one commit "goldens: regenerate for format v2". Both configurations green; TSAN job green.

### Wave 5 (cleanup and finish)

#### WP5.2 CRAM. Sonnet, M [R] (first package dropped if the schedule slips)
- **Owns:** `src/htsinterface/htsglue.c`, `thirdparty/htslib.cmake`, `packaging/bioconda/meta.yaml`, `.github/workflows/cmake.yml`, `tests/reader/CMakeLists.txt`, `tests/data/cram/` (new, small).
- **Cost, stated honestly:** 2 build dependencies (bz2, lzma) on every platform, CI changes on two OSes, recipe deps re-added, one CRAM fixture plus reference, 3 tests, and a support surface for reference errors.
- **Do:** enable bz2 and lzma in the bundled htslib (drop `--disable-bz2 --disable-lzma`); add `bzip2`, `xz` host deps back to the recipe; CI installs `libbz2-dev liblzma-dev` (ubuntu) / `xz` (macOS); set `CRAM_OPT_DECODE_MD`; reference from `REF_PATH`/`REF_CACHE` only. Missing reference: exit 1 with one line naming both variables; no network access.
- **Tests (label `reader`):** `cram:with_ref` (counts equal the same data as BAM), `cram:no_ref` (exit 1, message, `REF_PATH`/`REF_CACHE` unset, under a timeout), `cram:lzma_codec`.

#### WP5.3 stderr line and TTY progress. Sonnet, S
- **Owns:** `src/samstat.c`, `tests/cli/CMakeLists.txt`.
- **Do:** `standard-plots.md` section 6.
- **Tests (label `cli`):** `stderr_summary_line` (format regex, FASTQ variant without `% mapped`, `-p` suffix), `no_progress_without_tty`.

#### WP5.4 Docs and changelog. Sonnet, S
- **Owns:** `CHANGELOG.md`, `Readme.org`, `ChangeLog`, `NEWS`.
- **Do:** 2.4.0 entry from the collected PR lines; breaking changes listed first (JSON v2, protein FASTA unsupported, unknown input exit 1 per file, `-p` secondary/supplementary semantics, length now full read length, deprecated options); README usage and one screenshot-free description of the six plots.

#### WP5.1 Dead-code cleanup. Sonnet, S (after 5.3)
- **Owns:** `src/collect/collect.{c,h}`, `src/samstat.c`, `src/param/param.{c,h}`, `src/plot/plot.{c,h}`, `src/plot/plot_group.*`, `src/plot/plot_interval.*`, `src/plot/plot_test.c`, `src/plot/CMakeLists.txt`.
- **Do:** delete the legacy accumulators no longer read (Q x position matrix, clip-trimmed length, `_end` matrices, `len_norm`), `stat_collection_config_additional_plot` and its two calls (`samstat.c:171`, `:351`), the `collect_end`/`norm_len_plot` fields (`param.h:28-29`), the Plotly emitters, `plot_group.*`, `plot_interval.*`, `plot_test.c`, and `struct plot_data` once nothing uses it. Output unchanged.
- **Acceptance:** goldens byte-identical; `grep -r plot_data src/` is empty; `collect.c` is shorter.

#### Gate G5 (release, maintainer)
- `ctest` shows zero disabled tests (a still-disabled test means its feature was cut: delete it).
- Offline tarball build (as in WP0.1) for 2.4.0.
- Manual checks: one real Dorado BAM and one real Revio BAM (`decisions.md` E4); count how many of the real BAMs at hand lack MD (input to decision B2); report opened in Chrome and Firefox and printed to PDF once; `bench.sh` numbers on x86-64 Linux pasted into `docs/benchmark.md`.
- Maintainer tags `v2.4.0`; the bioconda autobump bot opens the version PR.

## 5. After 2.4.0 (not part of this workflow)

bio.tools entry; nf-core module (template `modules/nf-core/cramino`, after bioconda 2.4.0); MultiQC module-request issue; calibrated checks and caption ranges (2.5).
