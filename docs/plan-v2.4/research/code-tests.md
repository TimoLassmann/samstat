# Test infrastructure review and plan for v2.4 tests

Scope: tests/, .github/workflows, packaging/ at branch planning/v2.4 (2.3.0 merged). Read-only research; nothing modified.

## 1. How the suite works today

- Pure CMake, no external test framework. `tests/CMakeLists.txt` globs `dev/*.sam|*.fa` into `smoke:<file>` tests, then auto-`add_subdirectory()`s every `tests/<area>/` that has a CMakeLists.txt. A new area = new directory, nothing else to register. Shared variables: `SAMSTAT_TEST_DATA_DIR` (dev/), `SAMSTAT_TEST_OUT_DIR` (build/tests/out), `SAMSTAT_FIXTURE_DIR` (tests/data).
- Areas / labels: `build`, `cli`, `collect`, `fastx`, `golden`, `json`, `multiqc`, `packaging`, `reader`, `report`, plus `smoke`. Label = area name, set per test with `set_tests_properties(... LABELS x)`.
- Helpers in `tests/cmake/*.cmake`, run via `cmake -P` with `-D` args (portable, no Python/jq needed):
  - `run_smoke` (exit 0, non-empty HTML, `</html>` in last 64 bytes), `run_json` (+`CHECKS` or `GOLDEN`), `run_collect` (extra args, JSON checks, HTML regex present/absent), `run_threads` (-t1 vs -t4 identical), `run_golden`/`compare_golden`, `run_expect(_error)`, `run_write_failure`, `run_mqc`, `run_multiqc_tool`, `check_html`, `check_schema_keys`, `check_packaging`, `check_tld_pin`.
  - `check_json.cmake`: dotted paths, `path=value`, `path#len=N`, bare `path` = exists. Values cannot contain `;` or `=`. Numbers compared as CMake normalises them (50.00 == 50.0).
- Declarative wrappers: `collect_test()`, `json_value_test()`, `json_golden_test()`, `reader_json_test()` take `name input CHECKS...`. Hand-computed expectations live next to the fixture description in `tests/data/README.md` (expected values derived from the SAM spec, not from samstat output - keep this rule).
- One C unit test: `tests/reader/unit_parse_alignment.c`, linked against htsglue.c, sam_bam_parse.c, libhts and tld; compiled with ASAN/UBSan when `SAMSTAT_SANITIZE=ON`.
- Goldens (`tests/golden/expected/<input>.json`, 428 KB, ~one per dev/ and tests/data input) compare full `--json` output minus `meta`. `update-goldens` target regenerates them. `run_golden` silently skips inputs whose report exceeds `MAX_BYTES` (4 MB) - i.e. long reads are explicitly NOT golden-tested. `compare_golden` only reports the first differing top-level key (no numeric tolerance; exact string(JSON EQUAL)). A second, older golden set lives in `tests/json/golden/` updated by a manual cmake command (two mechanisms; merge them in 2.4).
- Fixtures: `tests/data` 140 KB (tiny hand-written SAM, small.bam, gz FASTA/FASTQ, truncated/garbage inputs). `dev/` is 4.1 MB of htslib-derived SAMs; `ce#large_seq.sam` (2.1 MB), `xx#large_aux.sam` (905 KB) and `ce.fa` (1.1 MB) are the only big ones. Note dev/ is the htslib test corpus, not curated for samstat; smoke/json/collect/golden all glob it, so adding a file to dev/ silently adds 4+ tests (and requires a golden).
- Large-report guard: `run_json`/`run_collect` skip JSON parsing above 4 MB and check only framing, because `string(JSON)` is slow. Any v2.4 JSON redesign that makes reports small (long-read binning) makes this guard unnecessary and brings long-read inputs under full value checks.
- Notable gaps: no 100k+ read test, no performance measurement, no memory check, no long-read data (max read 1 kb in `read_1024.sam`, `long_after_short.sam` 10 kb file), no MM/ML/supplementary-chain test beyond `supplementary.sam`, no stdin test, no offline-report test (HTML references the Plotly CDN), no HTML rendering test (only string regexes).

## 2. CI

`cmake.yml`: ubuntu-release, macos-release, ubuntu-asan (`-DSAMSTAT_SANITIZE=ON`, leak detection on), ubuntu-multiqc (pip multiqc==1.25.1, `-DSAMSTAT_REQUIRE_MULTIQC=ON`, `--schedule-random --repeat until-fail:2` to catch shared output dirs). All run `ctest -j4 --no-tests=error`. `release.yml`: on tag `v*`, build+ctest, make source tarball (with tld submodule) via `packaging/make_source_tarball.sh`, upload sha256. Packaging: draft `packaging/bioconda/{meta.yaml,build.sh}` (placeholder sha256, test = `samstat -v`), `packaging:metadata` and `packaging:tld_pin` CTest checks keep versions in sync.
Rules to preserve: random-order safe, each test writes to its own `OUT_DIR/<name>`, no network in tests, no tool beyond CMake + samstat (+optional multiqc).

## 3. Recommendations for the v2.4 implementation workflow

### 3.1 General principles (simplicity)
1. Keep the cmake-script style. Do not introduce pytest/Catch/etc. Python only for fixture generation (run as `uv run python`), never in the test run.
2. Add areas by adding directories: `tests/longread/`, `tests/perf/`. Reuse `collect_test()`-style wrappers; move that function into a shared `tests/cmake/helpers.cmake` included from `tests/CMakeLists.txt` if two areas need it (small refactor, justified by reuse).
3. Hand-computed expectations only for small fixtures; generated fixtures get expectations from the generator's own ground truth (it writes a sidecar `.truth.json` or the README table), never from samstat output.
4. Every v2.4 feature that changes a plot/JSON key gets exactly one value test and one golden; no per-feature test explosion.

### 3.2 Long-read fixtures (commit generator + outputs, small)
Add `tests/data/longread/` with a generator `tests/tools/make_longread_fixtures.py` (stdlib only, seeded `random.Random(24)`, run with `uv run python`; the outputs are committed so CI never runs it). Output budget: total under ~300 KB, each file under 100 KB.

| Fixture | Content | What it pins |
|---|---|---|
| `ont_like.sam` | ~40 reads, lengths log-normal 500 bp - 30 kb on a 60 kb synthetic reference, error profile ~8% (mostly indel, homopolymer deletions), MAPQ 0-60 mix, some soft-clipped ends, `MD`+`NM`, qualities Q5-Q25 | N50, log-length histogram, identity, error profile by MAPQ at long reads |
| `hifi_like.sam` | ~30 reads, 10-20 kb, ~0.5% error, Q20-Q60 (HiFi quals reach Q60+; tests the Q41 cap decision), `MM`/`ML` tags with 5mC on CpG (tag must parse and be ignored or counted, no crash) | high-accuracy regime, MM/ML tolerance |
| `supp_chain.sam` | one 20 kb read split into primary + 2 supplementary (flag 2048, `SA:Z` tags, hard clips on supplementary) + 1 secondary (256) | secondary/supplementary counting, read length from hard-clip + aligned (must not double count length or bases) |
| `ultralong.sam` | 2-3 reads of 1 Mb, SEQ `*`-free but low complexity so it compresses; ~3 MB raw. Keep OUT of git if >200 KB: generate at test time with a cmake-driven `uv run python` ONLY when env `SAMSTAT_PERF=1` | the 1 Mb -> ~100 MB JSON regression: assert JSON/HTML size is below a bound (e.g. < 2 MB) and independent of read length |
| `ont_like.bam` | `samtools view -b` of ont_like.sam (htslib 1.16 compatible; commit binary, ~15 KB) | BAM path with long CIGARs (> 65535 ops: put one read with a CG-tag long CIGAR, i.e. CIGAR placeholder `kSmN` plus `CG:B:I` tag) |
| `long.fq` / `long.fq.gz` | 20 reads 1-50 kb, FASTQ only (no alignment), plus one unterminated last record | FASTQ long-read path, N50, no-alignment mode, tld unterminated-record workaround |
| `unsorted_long_unmapped.sam` | unaligned (flag 4) uBAM-style long reads with `*` ref, `MM/ML` tags | unaligned long-read BAM input (very common: Dorado/ONT uBAM) |

Tool availability: `samtools` and `minimap2` exist on the maintainer's Mac but not as a guaranteed dependency; the generator must therefore write SAM directly (compute CIGAR/MD from its own simulated edits) and only the one-off `.bam` conversion uses samtools or `htslib test_view` (the existing README says BAMs were made with htslib's test_view - keep that, it avoids a new tool).

Assertions per long-read test (via `collect_test`/`json_value_test`): total_reads, mapped, supplementary/secondary counts, total bases, N50 (hand-computable from the generator's lengths; the generator prints them into `tests/data/longread/README.md`), length-histogram bin count is O(log length) not O(length), error-profile arrays have a bounded number of bins (e.g. `#len<=100`), mean identity within +-0.5% of the generator truth (use a range check; extend `check_json.cmake` with `path>=v` and `path<=v` operators - a ~10-line change and the one test-infra feature really needed since long-read statistics are not exact integers). Also `HTML_ABSENT "Infinity|NaN"` for every long-read run, and `HTML_PRESENT` for each of the standard plot titles (this enforces the "same plots every time" principle: one test asserts the fixed list of plot headings appears in every report, short and long).

### 3.3 Fixture size and repo hygiene
- Hard rule: new tracked fixtures total <= 500 KB; nothing > 100 KB except via generator. Document in `tests/data/longread/README.md` (same format as the existing README: purpose, expected values, how derived).
- Do NOT add long-read files to `dev/` (it is globbed into 4 test families and 2 golden sets).
- Prefer a bigger input generated on demand (see perf) over a big committed file. Do not use Git LFS.
- macOS is case-insensitive: no two fixtures differing only by case (existing note).

### 3.4 Performance tests: a `perf` label, off by default
Goal: catch regressions the maintainer cares about (speed, memory, report size) without flaky timing gates.
- New `tests/perf/CMakeLists.txt`. Tests are registered only when `-DSAMSTAT_PERF_TESTS=ON` (default OFF) so normal `ctest` and CI stay fast; labels `perf`. Run with `ctest -L perf`.
- Data: `tests/cmake/gen_perf_data.cmake` -> `uv run python tests/tools/make_perf_data.py --reads N --kind short|long --out <build>/perf/x.sam` (seeded, streaming writer; 1M x 150 bp short ~ 600 MB is too big for CI, so default 200k short reads and 2k long reads x 20 kb, scalable by `-DSAMSTAT_PERF_READS=`). A fixture test (`FIXTURES_SETUP perf_data`) creates the files once in the build tree; perf tests declare `FIXTURES_REQUIRED`. Fallback if `uv` absent: the test is skipped via `SKIP_RETURN_CODE 77`.
- Measurements: wrap in `tests/cmake/run_perf.cmake`: take timestamps with `string(TIMESTAMP ... %s)` has only second resolution, so use `cmake -E time`-free approach: `execute_process(COMMAND samstat ...)` between `cmake_host_system_information`-free timers is imprecise; better, let samstat itself print wall time and peak RSS with a `--verbose` line (it already has `--verbose`), or use `/usr/bin/time -l` (macOS) / `-v` (Linux). Simplest robust choice: one tiny shell wrapper `tests/perf/measure.sh` (POSIX sh, uses `/usr/bin/time`) invoked from cmake on UNIX. UNVERIFIED which of these samstat prints today.
- Gates (generous, to avoid flakiness; compare ratios not absolutes):
  1. Report size: long-read JSON and HTML < 2 MB regardless of read length (deterministic, not flaky - this one can run in the normal suite using `ultralong.sam`-style small data).
  2. Scaling: time(2N reads) / time(N reads) < 2.6 (catches quadratic behaviour such as the tld_strbuf_resize problem).
  3. Thread scaling (if v2.4 adds intra-file parallelism): `-t 4` wall time < 0.75 x `-t 1` on a single large BAM, only asserted when `nproc >= 4`, otherwise skipped (77). Also the existing `run_threads` equality check must be extended to a single multi-chunk BAM: results at `-t1` and `-t4` identical (this one is deterministic and belongs in the normal `cli`/`reader` label with a ~2 MB generated BAM, since chunking bugs are correctness bugs).
  4. Peak RSS < a fixed ceiling (e.g. 500 MB for the 2k x 20 kb input); streaming/`-p` sampling must keep memory flat.
  5. Absolute throughput is only RECORDED (printed, written to `build/perf/results.json`), not gated.
- CI: a separate `perf` job, `workflow_dispatch` + nightly cron on ubuntu-latest, Release build, uploads `results.json` as an artifact. Not required for PRs. Do not benchmark on macOS runners (noisy).
- Comparison with other tools (FastQC, NanoPlot, samtools stats) is a manual benchmark script for the paper/README, not a test.

### 3.5 Golden updates
- Expect the v2.4 JSON format change (new keys: length histogram, N50, identity; fewer/rebinned position plots) to invalidate every golden. Procedure for the implementation workflow:
  1. Land the format change; run `cmake --build build --target update-goldens`; commit goldens in a separate commit "goldens: regenerate for format v2" so the diff is reviewable; bump `format_version` in JSON, `docs/json-format.md`, `docs/samstat.schema.json`; `check_schema_keys` already ties schema to output.
  2. Review rule: a golden diff in a commit that is not labelled as a format change must be rejected.
  3. Merge the two golden mechanisms (`tests/golden/expected` and `tests/json/golden`) into the first; make `update-goldens` the single path.
  4. Lift the 4 MB `MAX_BYTES` skip once reports are small; add the long-read fixtures to the golden candidate glob (`tests/data/longread/*.sam|bam|fq`; extend the glob in `tests/golden/CMakeLists.txt` and handle FASTQ).
  5. Improve `compare_golden`: print the first differing path (not just top-level key) and allow a relative tolerance (1e-6) for floating point members (identity, mean quality), since long-read stats are floats and platform (macOS vs ubuntu) libm differences would otherwise cause CI-only failures. UNVERIFIED that this occurs today; current goldens are mostly integer counts.
  6. A "dev/ goldens are too many" check: consider restricting goldens to tests/data and a handful of dev/ files to cut the maintenance cost (428 KB of JSON diffs on each format change).
- HTML is not golden-compared and should stay that way (Plotly JSON embedded; brittle). Keep regex checks for presence of the fixed plot titles and absence of `NaN|Infinity|undefined`.

### 3.6 Other tests to add as features land (one each, minimal)
- Offline report (if Plotly is bundled/inlined): `HTML_ABSENT "cdn.plot.ly|https?://"` for the default report; size bound check.
- stdin streaming (if added): `samtools view -h x.bam | samstat -` equivalent via `run_expect` with `INPUT_FILE`; output name convention test; exit codes 0/1/2 unchanged.
- Native MultiQC module (if built): reuse `SAMSTAT_REQUIRE_MULTIQC` CI job, add a case with the long-read fixture and a short-read one in the same MultiQC run.
- Exit codes: keep existing `run_expect_error` tests; add truncated long BAM (EOF block missing) -> exit 2, no partial report claimed as success.
- Sanitizer: all new long-read tests carry no label exclusions, so the ASAN job exercises them (large CIGAR, 1 Mb reads are exactly where overflows hide; `read_1024.sam` shows this class of bug existed).
- Fuzz-lite: a loop test that truncates `ont_like.bam` at 5 offsets (cmake `file(READ ... LIMIT)`/`file(WRITE)`) and asserts exit code in {1,2}, never a crash (signal). Cheap, no fuzzer dependency.

### 3.7 CI changes
- Keep the 4 existing jobs. Add: (a) `perf` job (cron + manual) as above; (b) nothing else. Optionally pin `actions/checkout@v4` as is. Add ctest `-LE perf` is unnecessary because perf tests are not registered by default.
- Windows: not supported today; do not add.
- Release workflow: add the packaging checks already in CTest; after tag, the bioconda sha256 must be filled manually (draft recipe); a CTest check could fail if the placeholder sha256 survives while `PROJECT_VERSION` has a matching git tag - low value, skip.
- Update `packaging/bioconda/meta.yaml` version and test commands (`samstat -v`, add a `samstat --help` and a run on a tiny bundled file is not possible in bioconda `test:` without files; `source_files: tests/data/small.sam` is possible - UNVERIFIED, check conda-build docs https://docs.conda.io/projects/conda-build/en/stable/resources/define-metadata.html#test-section before using).

### 3.8 Suggested task split for the implementation workflow
1. Infra: `check_json` range operators; shared helpers file; unify golden mechanisms; improve compare_golden diagnostics. (small, do first, others depend on it)
2. Fixture generator + committed long-read fixtures + README with hand/generator-derived truths.
3. Long-read value + golden tests (needs the new JSON keys; written against the spec in docs/json-format.md before the implementation lands, so they fail first).
4. Perf label, data generator, `perf` CI job.
5. Golden regeneration commit after each format-affecting feature; one agent owns it to avoid merge conflicts in tests/golden/expected.
Conflict note: `tests/CMakeLists.txt` and `tests/golden/CMakeLists.txt` are shared edit points; assign them to a single agent. New areas (`tests/longread`, `tests/perf`) are conflict-free because of the directory auto-glob.

## 4. Unverified / external
- External tool features were not looked up; nothing in this note depends on a cited URL except the conda-build `test:` section (flagged UNVERIFIED).
- Whether samstat currently prints elapsed time / peak RSS under `--verbose` was not checked (UNVERIFIED).
- The current CTest count of 314 was taken from the task description, not re-run.
