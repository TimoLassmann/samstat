# SAMStat tests

The suite is plain CTest: every test runs the `samstat` binary of the build tree
through a small CMake script (`cmake -P`, helpers in `tests/cmake/`) or is a unit
test program built from `tests/<area>/`. It needs nothing beyond CMake, a C
compiler and, for the optional `multiqc:run` test, an installed `multiqc`.

## Running the suite safely

Configure and build in a directory of its own, outside backed-up folders if you
can (a scratch disk, `/tmp`), then run CTest from it:

```sh
cmake -S . -B ../samstat-build
cmake --build ../samstat-build -j4
(cd ../samstat-build && ctest --output-on-failure -j4)
```

- Use `-j1` to `-j4`. The tests are small, but several of them start samstat
  with up to 8 threads; `-j16` on a laptop means well over 100 threads at once
  and gains little.
- `ctest -L <label>` runs one area (`cli`, `reader`, `multiqc`, `report`, ...);
  `ctest -R <regex>` runs single tests.
- Everything a test writes goes below `<build>/tests/out/` (and generated inputs
  below `<build>/tests/<area>/`). Removing the build directory removes it all.
- No test uses the network. The CRAM tests use the references in `tests/data/cram`
  through `REF_PATH`; `cram:no_ref` (no reference at all) is registered only for
  the bundled htslib, which is built without libcurl. `multiqc:run` runs MultiQC
  with its version check off and with `HOME` and `XDG_CONFIG_HOME` inside the
  build tree.

## Rules for every test

These rules exist because a careless test once took down the maintainer's
laptop (`-t 100000`) and another broke Time Machine backups (a directory path
of about 1000 characters in the build tree). Review new tests against them.

1. **Small inputs.** Committed fixtures live in `tests/data` and the whole
   directory stays within 1 MB (`packaging:fixture_budget` checks it). Larger
   inputs are generated at configure or test time into the build tree and stay
   at a few MB at most (the 1 Mb single-read tests are the upper end). Lower a
   built-in limit with the test-only environment variables (for example
   `SAMSTAT_MAX_READ_LEN`, `SAMSTAT_BATCH_RECORDS`) instead of making an input
   big enough to reach it.
2. **Few threads.** Run samstat with `-t 8` at most. Tests of the option's own
   limit (`cli:threads_max_ok`, `-t 64` on a 20-read file) are the only
   exception; never test with a large value that would be accepted
   (`-t 100000`): if the cap broke, the test itself would start that many
   threads.
3. **Timeouts.** Every `execute_process` that runs samstat (or MultiQC) has a
   `TIMEOUT` (120 s for samstat), so a hang cannot block the run; a test that
   runs longer than CTest's default allows sets the CTest `TIMEOUT` property.
4. **Build tree only.** Write only below `SAMSTAT_TEST_OUT_DIR`
   (`<build>/tests/out`) or `CMAKE_CURRENT_BINARY_DIR`. Pass `-d @OUTDIR@`
   (or `-d ${D}`) to every samstat run, including runs that are expected to
   fail, so that a regression never writes a report next to a source-tree
   input. The helpers in `tests/cmake` refuse an `OUTDIR` that is not below
   `<build>/tests/out` (`outdir_guard.cmake`).
5. **No deep or long paths.** Directory names are short and a few levels deep;
   `run_expect.cmake` rejects a `MKDIR` entry over 100 characters or 4 levels.
   Path length limits are tested in unit tests, not by creating the path.
   The last test, `zz_path_check`, fails if anything below `tests/out` is
   longer than 150 characters or deeper than 8 levels.
6. **No network.** Tests that could reach the network with some build
   configuration (system htslib with libcurl, MultiQC's version check) are
   registered only where they cannot, or run isolated.
7. **Clean up what is large or awkward.** Remove large generated files at the
   end of the test. A test that makes a directory read-only restores the
   permissions at its start as well as its end and has a `FIXTURES_CLEANUP`
   test that restores them, so a killed run cannot leave a directory that
   `rm -rf build` fails on.
8. **No test runs Python, the benchmark or anything outside CMake.** The
   long-read fixtures are generated once with
   `tests/tools/make_longread_fixtures.py` and committed.

## Benchmark

`tests/perf/bench.sh` is not a CTest. It needs `BENCH_DIR` set to a scratch
directory outside any git work tree, prints the expected disk use before it
generates anything and caps the input sizes and thread counts unless
`BENCH_FORCE=1` (see the header of the script and `docs/benchmark.md`).
