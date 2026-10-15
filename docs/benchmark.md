# Benchmark

`tests/perf/bench.sh` measures wall time, peak RSS, Mreads/s and Mbases/s. It is not a CTest and gates nothing: timing on shared machines is noisy. Results are pasted below.

## Running it

```sh
cmake -S . -B ../samstat-build -DCMAKE_BUILD_TYPE=Release && cmake --build ../samstat-build -j
BENCH_DIR=/tmp/samstat-bench sh tests/perf/bench.sh ../samstat-build/src/samstat [path/to/samstat-2.3.0]
```

Keep the build and `BENCH_DIR` outside the repository: the script refuses to run without `BENCH_DIR`, or with a `BENCH_DIR` inside a git work tree, and caps the input sizes below unless `BENCH_FORCE=1` (the published results used `BENCH_FORCE=1 BENCH_SHORT_READS=2000000 BENCH_ONT_READS=150000`). The script generates SAM inputs with `tests/tools/make_longread_fixtures.py --reads N` (run through `uv run python`):

| Input | Default size |
|---|---|
| short reads, 150 bp | 200 000 reads (`BENCH_SHORT_READS`, cap 200 000) |
| ONT-like reads, 500 bp to 30 kb | 10 000 reads, about 70 Mb (`BENCH_ONT_READS`) |
| 1 Mb reads | 8 reads (`BENCH_HUGE_READS`) |
| FASTQ.gz | the short reads, converted with awk and gzip |

Generation time is printed on its own line and is not part of the timed runs. Each input is run at `-t 1`, `-t 4` and `-t 8` (`BENCH_THREADS`) under `/usr/bin/time` (`-l` on macOS, `-v` on Linux). The header lists the machine, the samstat version, the htslib version and whether the bundled htslib was built with libdeflate. When a second binary is given (2.3.0), it is timed on the same files.

## Targets

Targets are ratios:

- `-t 4` is at least 2.0x and `-t 8` at least 3.5x the `-t 1` rate of the same build, for inputs of one file.
- FASTQ.gz is at least 1.5x the 2.3.0 rate.
- 8 reads of 1 Mb at `-t 8` stay below 1 GB peak RSS.

## Results

Fill in one block per run. Build: bundled or system htslib, libdeflate yes or no.

| Machine | OS | Build | Input | -t | Wall (s) | Peak RSS (MB) | Mreads/s | Mbases/s |
|---|---|---|---|---|---|---|---|---|
| Apple M3 Max | macOS | bundled htslib | short reads SAM | 1 | | | 0.83 | |
| Apple M3 Max | macOS | bundled htslib | short reads SAM | 4 | | | 1.54 | |
| Apple M3 Max | macOS | bundled htslib | short reads SAM | 8 | | | 2.86 | |
| Apple M3 Max | macOS | bundled htslib | ONT-like SAM | 1 | | | | 97 |
| Apple M3 Max | macOS | bundled htslib | ONT-like SAM | 4 | | | | 208 |
| Apple M3 Max | macOS | bundled htslib | ONT-like SAM | 8 | | 553 | | 416 |
| Apple M3 Max | macOS | bundled htslib | FASTQ.gz | 1 | | | 0.95 | |
| Apple M3 Max | macOS | bundled htslib | FASTQ.gz | 4 | | | 1.67 | |
| Apple M3 Max | macOS | bundled htslib | FASTQ.gz | 8 | | | 1.54 | |

The FASTQ.gz rate does not improve from `-t 4` to `-t 8`. Empty cells were not recorded.
