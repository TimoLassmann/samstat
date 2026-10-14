# Benchmark

`tests/perf/bench.sh` measures wall time, peak RSS, Mreads/s and Mbases/s. It is not a CTest and gates nothing: timing on shared machines is noisy. The maintainer runs it and pastes the results below.

## Running it

```sh
cmake -S . -B ../samstat-build -DCMAKE_BUILD_TYPE=Release && cmake --build ../samstat-build -j
BENCH_DIR=../samstat-bench sh tests/perf/bench.sh ../samstat-build/src/samstat [path/to/samstat-2.3.0]
```

Keep the build and `BENCH_DIR` outside the repository. The script generates SAM inputs with `tests/tools/make_longread_fixtures.py --reads N` (run through `uv run python`):

| Input | Default size |
|---|---|
| short reads, 150 bp | 2 000 000 reads (`BENCH_SHORT_READS`) |
| ONT-like reads, 500 bp to 30 kb | 150 000 reads, about 1 Gb (`BENCH_ONT_READS`) |
| 1 Mb reads | 8 reads (`BENCH_HUGE_READS`) |
| FASTQ.gz | the short reads, converted with awk and gzip |

Generation time is printed on its own line and is not part of the timed runs. Each input is run at `-t 1`, `-t 4` and `-t 8` (`BENCH_THREADS`) under `/usr/bin/time` (`-l` on macOS, `-v` on Linux). The header lists the machine, the samstat version, the htslib version and whether the bundled htslib was built with libdeflate. When a second binary is given (2.3.0), it is timed on the same files.

## Targets

Targets are ratios, measured on x86-64 Linux after WP3.1a (intra-file threading, WP4.1, for the thread scaling):

- `-t 4` is at least 2.0x and `-t 8` at least 3.5x the `-t 1` rate of the same build.
- FASTQ.gz is at least 1.5x the 2.3.0 rate.
- 8 reads of 1 Mb at `-t 8` stay below 1 GB peak RSS.

Until then 2.3 runs one thread per input file, so a single file does not scale with `-t`.

## Results

Fill in one block per run. Build: bundled or system htslib, libdeflate yes or no.

| Machine | OS | Build | Input | -t | Wall (s) | Peak RSS (MB) | Mreads/s | Mbases/s |
|---|---|---|---|---|---|---|---|---|
| | | | | | | | | |
