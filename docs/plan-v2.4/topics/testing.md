# Topic: testing (final design for 2.4)

Base: the minimal proposal. Grafted from the researcher proposal: the mean-Q
calibration fixture, the MAPQ-bin sum invariant, the "point count does not grow
with read length" assertion, tests-first with a way to keep the suite green, and
`docs/benchmark.md`. Cut from both: perf CTest label, ratio gates, perf CI job,
float tolerance in `compare_golden`, `helpers.cmake`, uBAM fixture.

## Principle

The test suite is the contract the parallel implementation agents work to. It must
mechanically protect four promises to the researcher:

1. **Same plots every time.**
2. **Correct numbers.** N50, identity, mean Q and the MAPQ split match the truth.
3. **Bounded output and memory, whatever the read length.**
4. **The same answer for any `-t`.**

Keep the existing pure-CMake harness (`cmake -P` helpers, one directory per area,
auto-registered by the glob in `tests/CMakeLists.txt:36`). Python runs only to
generate fixtures, through `uv run python`, and never during `ctest`.

## Codebase facts relied on (verified by reading the code, 2026-10-14)

| Fact | Where |
|---|---|
| Test areas are auto-registered from `tests/*/` | `tests/CMakeLists.txt:36-39` |
| `dev/` is `SAMSTAT_TEST_DATA_DIR` and is globbed into the smoke and golden sets | `tests/CMakeLists.txt:9,14-16`, `tests/golden/CMakeLists.txt:17-19` |
| Golden glob covers `dev/*.sam`, `dev/*.fa` and `tests/data/*.{sam,bam}`. It does not cover FASTQ or subdirectories. | `tests/golden/CMakeLists.txt:17-19` |
| 4 MB skips | `run_golden.cmake:9-10,28`, `run_json.cmake:51`, `run_collect.cmake:34`, `run_expect.cmake:121` (NaN check) |
| A second golden set exists: 4 files, 48 KB | `tests/json/golden/` |
| `tests/golden/` is 428 KB | `du` |
| `run_threads` compares `-t1` and `-t4` across 8 copies of a file (8 files). Intra-file parallelism is therefore untested today. | `tests/cmake/run_threads.cmake` |
| `compare_golden` uses an exact `string(JSON EQUAL)` and reports no path | `tests/cmake/compare_golden.cmake` |
| JSON today writes integers only (`PRIu64`), with no floats | `src/report/json_report.c:328` |
| A truncated BAM test exists, for one file only | `tests/reader/CMakeLists.txt:139` |
| The multiqc CI job already runs `--schedule-random --repeat until-fail:2` | `.github/workflows/cmake.yml:93` |
| The ASAN job runs the full `ctest` | `.github/workflows/cmake.yml:61-70` |

## Decisions

### T1. Infrastructure first (wave 0, one agent, about 0.5 day)

1. **`check_json.cmake` gets `>=` and `<=`.** Exact `=` stays the rule for integer
   counts. Example: `summary.n50>=14850;summary.n50<=15150`.
2. **`compare_golden` names the first differing dotted path.** At present it only
   reports that the files differ. It does a recursive walk with `string(JSON ... GET/LENGTH/MEMBER)`.
   - No float tolerance. Instead, the JSON contract (owned by the reporting topic)
     prints every float with a fixed format such as `%.4f`.
   - Every float is derived from integer counts. Transcendental functions are
     replaced by lookup tables, for example the Phred-to-error-probability table
     for mean Q.
   - This makes the output byte-identical across glibc and macOS libm without
     tolerance code. A tolerance is added only if a CI-only diff ever appears.
3. **One golden mechanism.**
   - Move `tests/json/golden/*` into `tests/golden/expected/` and delete `tests/json/golden/`.
   - `update-goldens` becomes the only way to regenerate goldens.
   - Extend the golden glob to `tests/data/longread/*.{sam,bam,fq}`.
4. **A fixture budget check (packaging label, about 5 lines).** It fails if
   `tests/data` grows past 1 MB, because a multi-agent run will otherwise commit
   large files.
5. **Disabled-until-ready tests.** Wave 0 writes the feature tests against the spec
   in `docs/json-format.md` and registers each with `DISABLED TRUE`.
   - CTest reports disabled tests as "Not Run (Disabled)" and does not fail the run,
     so the suite stays green and CI commands are unchanged. No custom `pending` label is needed.
   - A feature package removes `DISABLED` from its own tests. Its acceptance gate is
     that those tests pass.

The 4 MB skips are removed in the package that bounds the per-position arrays
(not in wave 0). Removing them earlier would make `dev/` long-read inputs produce
~100 MB goldens.

### T2. Long-read fixtures: five small committed files

Location: `tests/data/longread/`. The generator is
`tests/tools/make_longread_fixtures.py`: stdlib only, seeded, run once with
`uv run python`. Its outputs are committed, so CI never runs Python. BAMs are made
with htslib `test/test_view`, as the existing `tests/data/README.md` does.

| File | Content | Pins |
|---|---|---|
| `ont_like.sam` and `ont_like.bam` (same reads) | ~40 reads, 500 bp to 30 kb, ~8% error, a mix of MAPQ values, soft clips, MD and NM. The BAM adds one read with more than 65535 CIGAR ops (CG tag). | N50, the log-length histogram, identity, bounded error-profile bins, the BAM path, the long-CIGAR path |
| `supp_chain.sam` | 1 primary, 2 hard-clipped supplementary and 1 secondary alignment | Read count, length and bases are not double-counted |
| `hifi_like.sam` | 15 kb reads, Q up to 60, MM/ML tags | Q>41 is not saturated. MM/ML tags are tolerated. It also feeds the modbase block if that ships. |
| `long.fq` | 20 reads of 1 to 50 kb with an unterminated last record. Two reads have identical Phred means but different error-probability means. | The FASTQ long-read path. Mean Q is computed from error probabilities, not by averaging Phred. |

Rules:

- **Where expected values come from.** The generator writes expected values from its
  own ground truth into `tests/data/longread/README.md`, in the style of
  `tests/data/README.md`. They are never copied from samstat output.
- **Hand-checkable cases.** Add small hand-computable cases next to the statistical
  ones, for example an N50 on 5 reads and the supplementary counts. These catch a
  generator that shares samstat's bug.
- **Size budget.** At most 300 KB in total and at most 100 KB per file.
- **Nothing goes into `dev/`.**

### T3. Value tests per fixture (`tests/longread/`)

**Integer counts, exact:**

- reads
- bases
- primary, secondary and supplementary counts
- the MAPQ bins `gt30`, `lt30`, `lt10` and `unmapped`, which must sum to the number of primary reads

**Ranges:**

- N50 within 1%. This is what the planned 1% log length bins can resolve.
- Median identity within 0.5 percentage points of the truth.
- Mean Q within 0.1 of the truth on `long.fq`.
- Max quality on `hifi_like.sam` is at least 60.

### T4. Bound on ultralong output (in the normal suite)

`run_ultralong.cmake` writes a single 1 Mb read SAM into `OUT_DIR` at test time.
The SAM has `string(REPEAT)` SEQ and QUAL, `CIGAR 1000000M` and `MD:Z:1000000`. It
is about 2 MB, is not committed, and needs no Python. The test asserts:

- The JSON is under 1 MB and the HTML is under 2 MB.
- The position arrays have the same `LENGTH` as those for `tests/data/small.sam`,
  or at most the fixed bin cap. The number of points must not grow with read length.

It runs in the ASAN job as well, which is where overflows in 1 Mb buffers would show up.

Peak RSS is not checked in CTest. `ulimit -v` cannot be used under ASAN, and
`/usr/bin/time` differs between macOS and Linux. RSS is measured by `bench.sh` (T7).

### T5. The fixed plot set (one test, highest value for the HTML)

- One CMake list variable holds the fixed plot titles in order. A second, shorter
  list holds the subset for FASTQ and unaligned input.
- The test runs one short BAM, `ont_like.bam` and `long.fq`. It asserts that the
  titles are present in the expected order.
- It asserts that the text `NaN|Infinity|undefined` is absent.
- It asserts that `https?://` is absent, which makes the report offline. This
  check is enabled when the presentation package lands.
- Changing the plot set means editing one line, which makes plot creep visible in review.
- HTML is never golden-compared.

### T6. The same answer for any `-t`, within one file

- A new case in `run_threads.cmake` (the existing case across 8 files stays).
- At test time it builds one multi-chunk SAM in `OUT_DIR`: the header of
  `tests/data/paired_small.sam` plus its body repeated until it holds at least 5
  reader chunks (more than 50,000 records at the current chunk size).
- It runs `-t 1`, `-t 3` and `-t 8` and requires the three JSON outputs to be
  byte-identical. `-t 3` catches partitioning off-by-ones.
- The same comparison runs on `ont_like.bam`.
- Races: the multiqc job already repeats in random order (`cmake.yml:93`), and ASAN
  runs this test. No new job.

### T7. Performance: one script, no gates

- `tests/perf/bench.sh` is POSIX sh. It is not a CTest and has no CI job.
- It calls the T2 generator with `--reads N` to make data in the build directory:
  about 2 M short reads and about 1 Gb of ONT-like reads.
- It runs samstat under `/usr/bin/time` with a small macOS/Linux wrapper, at
  `-t 1`, `-t 4` and `-t 8`.
- It prints wall time, peak RSS, Mreads/s and Mbases/s.
- The maintainer pastes the results, with the machine labelled, into
  `docs/benchmark.md`. This backs the paper's speed claim and the targets in
  `research/competitor-speed.md`, which are recommendations and not verified facts.
- Why there are no gates: timing on shared runners is noisy. The regressions that
  matter (output size, point count) are covered deterministically by T4.

### T8. Failure behaviour

- **Truncation loop (reader label).** Cut `ont_like.bam` at 5 offsets with
  `file(READ ... LIMIT)` and `file(WRITE)`. Each cut must exit with 1 or 2, not by a
  signal, and must leave no report behind. About 20 lines of CMake and no fuzzer.
  This is worth having because the reader is being rewritten.
- **CLI checks.** These are written by the CLI package as its value tests, only if
  the corresponding CLI change ships:
  - an unwritable `-o` fails before processing starts;
  - duplicate basenames with `-d` give an error or distinct outputs;
  - piping a file into `samstat -` gives the same JSON as reading the file;
  - a `.cram` input with a reference runs; a `.cram` input without one gives a clear exit 2 message.
- Exit codes 0, 1 and 2 are unchanged.

### T9. Goldens

- Goldens cover `tests/data/**`, the `dev/` files already covered and
  `tests/data/longread/*`. Add no more `dev/` goldens.
- **Regeneration commit.** Goldens are regenerated only in a dedicated commit,
  "goldens: regenerate for format vN". That commit also bumps `format_version` and
  updates `docs/json-format.md` and `docs/samstat.schema.json`.
  `check_schema_keys` already ties the schema to the output.
- **Merge rule.** A golden diff in any other commit is rejected.
- **Batching.** Format changes are batched so goldens are regenerated at most twice:
  1. after the bounded-bins and long-read statistics change;
  2. after the report and JSON reshaping.
- **Behaviour-preserving packages** (htslib threads, reader switch to htslib for
  FASTQ, htslib bump) must leave the goldens byte-identical. This is the strongest
  available check that nothing changed silently. If the htslib bump changes a
  golden, that is reviewed on its own as a separate package.

## How the multi-agent workflow keeps tests green

1. **Wave 0 (one agent).** Does T1, T2 and T5, writes the T3, T4 and T6 tests as
   `DISABLED`, and adds the T8 truncation loop. It merges green.
2. **Owner of shared test files.** One agent owns `tests/CMakeLists.txt`,
   `tests/golden/CMakeLists.txt` and `tests/golden/expected/`. Feature agents add
   only to their own `tests/<area>/` directory, which is conflict-free because of the
   glob.
3. **Each feature package:**
   - enables its own disabled tests;
   - adds at most one value test of its own;
   - must not edit goldens.
4. **Merge gate for each package.** All of the following must hold:
   - `ctest -j4 --no-tests=error` is green in Release;
   - it is green with `-DSAMSTAT_SANITIZE=ON`;
   - the golden diff is empty, or the commit is the labelled regeneration commit
     made by the golden owner.
5. **Order.** Format-changing packages are serialised. Behaviour-preserving
   packages can run in parallel.
6. **End of 2.4.** No test is left `DISABLED`. A test still disabled means its
   feature was cut, and the test is deleted.

## Explicitly out

- A perf CTest label, a perf CI job, and timing or ratio gates.
- Float tolerance in goldens, unless it is ever proven necessary.
- pytest, Catch, Unity, or Python at test run time.
- Git LFS and large fixtures. `ultralong.sam` is generated at test time and not committed.
- Windows CI, any expansion of the matrix, fuzzers, network access in tests.
- HTML goldens, and browser or screenshot tests.
- These fixtures: an unaligned uBAM (FASTQ covers that path), a separate ONT/HiFi
  BAM pair, and a modbase fixture set beyond `hifi_like.sam`.
- Real-data tests. A real Dorado or Revio check stays manual (see the risks).

## Affected files

- **Modified:**
  - `tests/cmake/check_json.cmake`
  - `compare_golden.cmake`
  - `run_golden.cmake`, `run_json.cmake`, `run_collect.cmake`, `run_expect.cmake` (4 MB skips)
  - `run_threads.cmake`
  - `tests/golden/CMakeLists.txt`
  - `tests/json/CMakeLists.txt` (it points to the merged goldens)
  - `tests/reader/CMakeLists.txt` (truncation loop)
- **New:**
  - `tests/longread/CMakeLists.txt`
  - `tests/cmake/run_ultralong.cmake`
  - `tests/data/longread/*` and its `README.md`
  - `tests/tools/make_longread_fixtures.py`
  - `tests/perf/bench.sh`
  - `docs/benchmark.md`
- **Deleted:** `tests/json/golden/`.

## Effort

| Item | Days |
|---|---|
| T1 infrastructure | 0.5 |
| T2 generator and fixtures (the CIGAR, MD and NM ground truth is most of the work) | 1 |
| T3 to T6 tests | 0.5 |
| T7 bench script | 0.25 |
| T8 truncation loop | 0.25 |
| T9 per golden regeneration | 0.25 each |
| **Total** | **about 3** |

## Risks

- **The generator could share a bug with samstat.** Mitigation: hand-computable cases
  (T2), and a one-off manual cross-check of counts against `samtools stats`.
  `samtools stats` defines identity differently, so only the counts are compared.
- **Synthetic data is not real data.** The tests pin counting, not realistic ranges.
  ONT `qs`, HiFi QUAL binning and typical identity ranges are UNVERIFIED (see
  `research/longread-metrics.md`). Before release, run one real Dorado BAM and one
  real Revio BAM by hand.
- **No perf gate.** A throughput regression can ship. This is accepted; `bench.sh`
  makes it visible on demand.
- **Removing the 4 MB skips slows the suite if bounding fails.** In that case T4
  fails first, which is intended.
- **htslib bump for MN.** If MM/ML is decoded, the htslib bump needed for MN
  (1.18 per `research/modbases.md`, not checked against the 1.16 headers) is
  UNVERIFIED. The bump is its own package, with its golden diff reviewed separately.

## DECISION NEEDED

1. **Float format in the JSON.** Option one is a fixed `%.4f`, computed from integers
   and lookup tables, with no tolerance code (recommended). Option two is a 1e-6
   tolerance in `compare_golden`.
2. **Peak RSS.** Option one is manual only, through `bench.sh` (recommended). Option
   two is a Linux-only CTest using `/usr/bin/time -v` with a ceiling such as 2 GB.
3. **Which `dev/` goldens to keep.** Keep all current ones (428 KB, regenerated twice
   in 2.4), or cut them to a short named list.
4. **Real-data check before release.** Is a manual run on one Dorado BAM and one
   Revio BAM required before tagging 2.4, and which files would be used? Licensing
   rules them out of CI.
5. **Ultralong ceilings.** JSON under 1 MB and HTML under 2 MB are proposed. They
   depend on the presentation decision between inline SVG and embedded Plotly:
   embedded Plotly adds 1.2 to 4.8 MB.
