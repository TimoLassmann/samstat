# SAMStat 2.4 topic design: cli-ux (final)

Base: the **minimal** proposal. Grafted from the researcher proposal: the layout is decided at the end of the run (which makes stdin work), the `-p` partial note goes in the report, and a test that checks `-h` against the README. Everything else in either proposal that does not serve `samstat file` doing the right thing has been cut.

## Guiding rule

- `samstat file` must work with no options.
- 2.4 has fewer options than 2.3.
- The output file names and the meaning of exit codes 0/1/2 are frozen. MultiQC globs, the goldens and the bioconda recipe depend on them.

## Verified state of 2.3.0

I checked these by reading the code and running `build-rel/src/samstat`. They correct some claims made in the proposals.

| Claim | Reality in 2.3.0 |
|---|---|
| `-l` and `--seed` are dead | `-l` is parsed and range-checked (`param.c:103`, `:225-231`), but nothing reads it outside `param.c`. **`--seed` is worse than dead.** It has an entry in the option table (`param.c:47`) but no `case` in the switch, so it falls through to `default: abort()`. **`samstat --seed 1 x.sam` crashes with rc 134.** That is a 2.3 bug the deprecation fixes. |
| An unknown extension aborts the run | Yes. `sanity_check_param` (`param.c:249-256`) makes it a usage error: **exit 2, and no file is processed**. Tests `cli:ext_unknown_only` and `cli:ext_unknown_with_good` (`tests/cli/CMakeLists.txt:92-96`) lock this behaviour in, so they must change. |
| `-v` output | It prints `SAMStat 2.3.0` (`PACKAGE_NAME PACKAGE_VERSION`). The bioconda test only runs `samstat -v` (`packaging/bioconda/meta.yaml:34`). |
| `-t` | The long name is `--nthreads`, not `--threads`. The help says "will only be used when multiple input files are present". |
| The output directory must exist | `param.c:200-201` errors if it does not. |
| The output path is checked late | It is built in `stat_report.c:261-272` after all records have been read. |
| Duplicate inputs | Only identical input *strings* are rejected (`param.c:260-266`). `run1/a.bam` and `run2/a.bam` with `-d` overwrite each other silently. |
| `-p` partial note | It exists on stderr (`samstat.c:230-231`) and in the JSON (`is_partial_report`, `json_report.c:99`). It is **absent from the HTML** (`stat_report.c` never reads the flag). |
| Parser | `getopt_long_only`, so `-dir`, `-peek` and `-norm` work as single-dash long options. Unambiguous prefixes are also accepted. |

## Decisions

### 1. Option set: 8 flags, none required

| Option | Status | Notes |
|---|---|---|
| `-o, --outdir DIR` | New primary name | Created if missing (`mkdir -p`). `-d`, `-dir` and `--dir` stay as silent aliases forever. |
| `-p, --peek N` | Keep | Help says "first N reads (not random; biased on sorted files)". |
| `-t, --threads N` | Keep. Add `--threads`, keep `--nthreads` as an alias | Means **total threads**. Default min(4, ncores). The single-file speedup is implemented by the perf topic. If that slips, the help text keeps the "multiple input files only" wording. |
| `--json` | Keep | |
| `--multiqc` | Keep | |
| `--verbose` | Keep | Debug only. |
| `-h`, `-v` | Keep | `-v` output stays **`SAMStat X.Y.Z`** (unchanged). Parsers use `sed 's/SAMStat //'`. Changing the string brings no benefit and breaks anyone already parsing it. |

### 2. Deprecations (accept, warn once, ignore; exit code unaffected)

| Option | 2.4 | Removal |
|---|---|---|
| `-l/-len N` | Value accepted without range checks. Prints `samstat: WARNING: -l is deprecated and ignored`. | 3.0 |
| `--seed N` | Same. **Fixes the 2.3 abort.** | 3.0 |
| `--norm` | Same. The layout is automatic (report/collect topics). Release note: a 2.3 user of `--norm` on short reads now sees the default layout. | 3.0 |
| `--plotend` | Depends on the report topic. If the end-anchored plots are dropped, it gets the same treatment. If they are kept, they become default-on for short reads. **No flag either way.** | 3.0 |

No `--mode short|long` override. If auto-detection misfires on real data, fix the heuristic. Add the override only after a real report of a misdetection.

### 3. Output location and naming

- **Names are unchanged:** `<input>.samstat.html`, `.samstat.json` and `.samstat_*_mqc.json`.
- **The default stays beside the input.** A CWD default is a 3.0 question.
- **Up-front checks, under 1 s and before any record is read.** Resolve every output path, then:
  - Check that each output directory is writable, with `access(dir, W_OK)`. Do not create a probe file, because that leaves debris on failure.
  - If it is not writable, fail **that input** with exit 1. The message names the path and suggests `-o`.
  - **No silent fallback to CWD.** Writing somewhere unexpected is worse than a clear error.
  - Output-path collision: if two inputs resolve to the same output, exit 2 before any processing, naming both inputs. This extends the existing duplicate check at `param.c:260`.
- **No** `--prefix`, `--name`, `--format`, and no always-on JSON.

### 4. Input detection, stdin, CRAM

- **The format comes from content, not the extension.** Use htslib `hts_open` + `hts_get_format`. This covers BAM, SAM, CRAM, FASTA, FASTQ, gzip and BGZF. The extension is used only as a tiebreaker for empty files.
  - This depends on the io topic replacing the tld FASTA/FASTQ reader and the `/tmp` copy in `prepare_fasta_fastq_input` (`samstat.c:257`) with htslib reading. code-io.md verified that htslib 1.16 reads FASTQ/FASTA, `.gz`, pipes and unterminated last records.
- **An unrecognised or unreadable input fails only that input** (exit 1). The other inputs still run, which matches how a missing file is already handled. Update `cli:ext_unknown_*` from RC 2 to RC 1, with the valid file's report present.
- **`-` means stdin.** It is allowed at most once per command line (a second `-` is exit 2).
  - The output is `stdin.samstat.html` (and siblings) in `-o`, or in the CWD if `-o` is not given. This keeps the `*.samstat*` glob pattern, and MultiQC shows the sample as `stdin`.
  - `<(...)` process substitution works too, because detection no longer needs an extension. Its output name is the substituted path's basename (e.g. `63`), so recommend `-` in the docs.
- **The layout (short/long) is decided at the end of the run, from the bounded length histogram, not at the start.** That needs no second pass, so stdin works. This requires the collector design (code-collect.md) to cooperate.
- **CRAM is accepted** (`.cram` or sniffed).
  - If a record has no MD tag, the mismatch/indel profiles are reported as **"not available: no MD tag"** in the report, with one stderr warning. They are never shown as silent zeros. How the 2.3 parser behaves without MD is UNVERIFIED; check this first, because it also affects long-read BAMs without MD.
  - **No `--reference` flag in the baseline.** htslib's `REF_PATH`/`REF_CACHE` and the `@SQ UR:` header cover it. If decoding fails, the error message names those variables. See the DECISION NEEDED item below.

### 5. Messages (stderr)

- **One final line per input, always printed:**
  `samstat: wrote x.bam.samstat.html (BAM, 12.3M reads, long-read layout: median 8.2 kb, 41 s)`
  - The layout choice is only known at the end, so it goes here and in the HTML header, not in a start line.
  - With `-p`, the line adds `, partial: first N reads`.
- **The HTML header gets the same layout line and the partial note.** The partial flag already exists, so this is a one-line addition in `stat_report.c`. It is shared with the report topic.
- **Progress:** a counter every ~1M reads, only when `isatty(STDERR_FILENO)`. Logs and CI stay clean.
- **No `-q`.** The output is one line per file, and users can redirect stderr.

### 6. Exit codes (unchanged meaning)

- 0 means OK.
- 1 means at least one input failed, while the others completed. In 2.4 this also covers an unrecognised input and an unwritable output directory.
- 2 means a usage error: a bad option or value, no inputs, an output collision, or `-` given twice.
- The codes are printed in `-h`.

### 7. `-h` text

```
samstat 2.4.0 - quality control for SAM/BAM/CRAM and FASTA/FASTQ files

Usage: samstat [options] <file|-> [file ...]

  -o, --outdir DIR   write reports to DIR (created if missing) [next to input]
  -p, --peek N       analyse only the first N reads (not random) [all]
  -t, --threads N    total threads [4]
      --json         also write <input>.samstat.json
      --multiqc      also write MultiQC files <input>.samstat_*_mqc.json
      --verbose      debug output
  -h, --help / -v, --version

Input format is detected from file content; '-' reads stdin (report: stdin.samstat.html).
Exit codes: 0 ok, 1 one or more inputs failed, 2 usage error.

Examples:
  samstat sample.bam
  samstat -o qc/ run/*.bam run/*.fastq.gz
  minimap2 -a ref.fa reads.fq | samstat -o qc/ -
```

Update `Readme.org` to match. It currently lists `--seed` and `--plotend` and omits `--norm` and `-l`. Add a deprecation note.

## Explicitly OUT

- `--mode short|long`.
- `--format`, `--prefix`, `--name`, `-q`, `--sub`, `--seed` as random sampling.
- A config file, shell completions, any interactive UI.
- A CWD fallback for unwritable directories.
- A change of the default output location.
- Always-on JSON.
- Configurable MAPQ bins.

## Affected files

- `src/param/param.c`, `param.h`:
  - option table (`-o`, `--outdir`, `--threads`)
  - deprecation cases, including the missing `OPT_SEED` case
  - removal of the `-l` range check
  - `mkdir -p` in place of the "does not exist" error
  - output-path resolution, writability and collision checks
  - remove the extension gate (`detect_file_type`, `sanity_check_param:233-256`)
  - new `print_help()`
- `src/samstat.c`:
  - per-file failure for unrecognised input
  - `-` handling
  - final summary line and TTY progress
  - the reader switch to htslib, shared with the perf and io topics
- `src/htsinterface/htsglue.c`: `hts_get_format`, stdin, CRAM.
- `src/report/stat_report.c`: layout and partial line in the header (report topic).
- `Readme.org`, `docs/`.
- `tests/cli/CMakeLists.txt`.

## Effort: about 3 days

1. Deprecations (including the `--seed` crash fix), `-o` with mkdir, `--threads`, the new `-h` and README: 0.5 day.
2. Up-front path resolution and the writability and collision checks: 0.5 day.
3. Content sniffing, `-`, CRAM, and per-file unknown input: 1.5 days. Partly shared with the io/perf reader rewrite, and dependent on it.
4. Final line, TTY progress, and the header lines: 0.5 day.

## Risks

- **Behaviour change.** An unknown input becomes exit 1 per file instead of exit 2 for the whole run. Note it in the release notes and update the two existing tests.
- **Sniffing edge cases.** Empty files, truncated gzip and `.bz2` must give clear per-file errors. Delegate detection to htslib and only handle the empty and missing cases ourselves.
- **stdin with `-p`.** Stopping early sends SIGPIPE upstream. That is normal, but it needs a test.
- **The always-on final line** may break any test that asserts empty stderr. Audit `tests/cli` and `tests/smoke`.
- **CRAM or BAM without MD.** Parser behaviour is UNVERIFIED. It must degrade to "not available", never to zeros.
- **`-t` honesty.** The help must not promise a single-file speedup until the perf work lands.

## Tests (label `cli`, in `tests/cli/CMakeLists.txt`)

1. `-o new/dir` is created and the report lands in it. `-d` and `-dir` give byte-identical output.
2. A read-only output dir gives exit 1 within 1 s, before reading. The message names the path. Other inputs still succeed.
3. `-o qc a/x.sam b/x.sam` gives exit 2 and names both inputs, with no files written.
4. `-l 100`, `--seed 1`, `--norm`: exit 0 (regression for the rc-134 abort), one warning line each, and output identical to a run without the flag (`--norm` is compared with the new default).
5. `cat x.bam | samstat -o D -` gives the same JSON as `samstat -o D x.bam` except for name and path fields. A BAM named `x.dat` and a gzipped FASTQ with no extension work. `- -` gives exit 2.
6. A batch of `small.sam` plus `unknown.xyz` gives exit 1, and `small.sam.samstat.html` exists. This replaces `ext_unknown_with_good`.
7. A `.cram` fixture with a local reference via `REF_PATH` works. Without MD, the report says "not available" and there is no crash.
8. `-v` matches `^SAMStat [0-9]+\.[0-9]+\.[0-9]+$`. `-h` contains the exit-code line. Every long option in `-h` appears in `Readme.org`.
9. The final line is present on stderr. The progress counter is absent when stderr is not a TTY. With `-p`, the HTML contains the partial note.
10. The existing goldens and the JSON and MultiQC filename tests pass unchanged, which guards the naming contract.

## DECISION NEEDED

1. **`--plotend`.** Deprecate it, or make its plots default-on for short reads? This follows the report topic's decision.
2. **Protein FASTA.** Dropping the tld reader in favour of htslib loses protein FASTA. Accept that loss with a clear per-file error? Shared with the io topic.
3. **CRAM reference.** Rely on `REF_PATH`/`REF_CACHE` only (recommended), or add `-r/--reference FILE`? perf-htslib.md notes that htslib 1.22+ no longer fetches references from EBI by default. That is UNVERIFIED here, and after the planned htslib bump it makes reference-less CRAM fail outright.
4. **stdin output name.** `stdin.samstat.html` is the recommendation. The alternative is a bare `samstat.html`, which breaks the `<name>.samstat` pattern.
5. **Default output location.** Keep beside the input in 2.4 (recommended), and leave CWD for 3.0.
