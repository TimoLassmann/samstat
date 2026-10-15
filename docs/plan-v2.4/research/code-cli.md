# CLI review and proposal for SAMStat 2.4

Scope: src/param/param.c, help text, Readme.org, docs/. Static read only (no build was available in this checkout, `samstat -h` output was taken from `print_help()` in param.c). Nothing external is cited; nothing UNVERIFIED except where marked.

## 1. Current CLI (2.3.0)

`samstat [options] file1 file2 ...` with getopt_long_only (so `-peek`, `-json` and `--json` all work).

| option | default | actual behaviour |
|---|---|---|
| `-p/-peek N` | unlimited | analyse exactly the first N reads (samstat.c:190, 373). Not random; first-N is biased for sorted BAM (only chr1 start). |
| `-l/-len N` | 250 | parsed and range-checked (1..10000, warns >500 "slow") but **no consumer**: `report_max_len` is never read (collect.c:594 is commented out). Dead option. Also NOT listed in help or README. |
| `-d/-dir` | next to input | output dir, must already exist (no mkdir). |
| `-t N` | 4 | threads; one file per thread, so no effect for a single file. |
| `--plotend` | off | extra read-end plots. |
| `--norm` | off | extra plots in 1% read-length bins. Listed in help, missing from README. |
| `--json` | off | extra output `<in>.samstat.json`. |
| `--multiqc` | off | extra output `<in>.samstat_*_mqc.json` (4 files). |
| `--seed N` | 0 | parsed, **never used** (`->seed` has no reader; subsampling `--sub` is commented out and sanity-checked against a hard-coded 1.0). Dead option. |
| `--verbose` | off | progress messages. |
| `-h`, `-v` | | help, version. |

Positional: one or more files; type from extension only (.sam .bam .fa .fasta .fna .fq .fastq + .gz for fasta/fastq). Unknown extension is a hard error for the whole run (before any file is processed). Missing file is reported at processing time (exit 1, others continue). Exit codes 0/1/2 (documented in README, not in `-h`).

Output: `<input>.samstat.html` beside the input (e.g. `reads.fastq.gz.samstat.html`), or `<dir>/<basename>.samstat.html`.

## 2. Problems, ranked by cost to a researcher

1. **No CRAM, no stdin, no extensionless input.** `.cram` is silently unsupported (error "could not recognize type"), although htslib reads it. `samtools view -h x.bam | samstat -` and process substitution (`<(...)`) cannot work because type is decided by the file name. Pipelines (minimap2 | samtools sort) cannot tee into QC. Fix: sniff magic bytes (BGZF 1f 8b + BAM\1 inside, `@` header / tab fields for SAM, `>` / `@` for FASTA/FASTQ, `CRAM` magic) and accept `-` for stdin; extension becomes a hint, not a gate. Also handles `.bam` files misnamed, `.fq.bz2`/`.fastq.zst` clearly rejected by content.
2. **Writes next to the input.** Fails (opaque fopen error after the whole file is processed) on read-only directories, e.g. shared data archives, which is the usual place for BAMs. Error comes at the END of a possibly hour-long run. Fix: probe writability of the output location up front, fail in <1 s; and default to the current working directory, not the input directory (FastQC writes next to the input by default and has the same complaint; the pipeline convention is `-o`).
3. **Dead options (`-l`, `--seed`) and a misleading one (`-t`).** `-t 4` default with help text "only used when multiple input files" invites users to pass `-t 16` for one BAM and get no speedup. Dead options are a trust problem and they were visibly removed from the help text, but still parse silently. Fix: remove `-l` and `--seed` (accept and ignore with a one-line deprecation warning for one release to keep scripts working). Make `-t` mean "threads" for real (htslib BGZF decompression threads via `hts_set_threads` for one file is the cheap 2-3x; the larger parallel-parse redesign belongs to the performance plan, not the CLI).
4. **Short/long reads are not distinguished by the CLI at all.** No option exists and none should: the report must adapt automatically (length binning, log-length histogram, N50). Today `--norm` (1% bins) is effectively the "long read" switch and is off by default, so a first-time nanopore user gets per-position plots that are meaningless beyond a few kb. Fix: auto-select from data (read-length distribution of the first N reads, e.g. median > 1 kb or max > 10 kb => long-read presentation), with the choice printed in the report header. No flag; at most one override (`--mode short|long`) as an escape hatch, only if the heuristic proves wrong in tests.
5. **Output naming.** `sample.bam.samstat.html` is unambiguous but ugly, and collides when `-d` is used with two inputs with the same basename in different directories (`run1/a.bam`, `run2/a.bam` -> second overwrites first silently). Fix: detect output-name collisions up front and error (cheap). Keep `<input>.samstat.html` name for backward compatibility (MultiQC/nf-core globs and goldens depend on it).
6. **Opt-in outputs.** `--json` and `--multiqc` are off by default, so a pipeline author must know two flags. Since JSON is the machine interface and MultiQC files are tiny, an alternative is "always write `.json`". Not recommended: it changes the file set users see. Recommended instead: a single `--json` stays; `--multiqc` is documented as the nf-core/MultiQC flag. Final answer depends on the plan for a native MultiQC module (see section 5).
7. **Unnecessary strictness.** One unrecognised extension aborts all files; a duplicated path is an error (fine); `-d` requires an existing directory (should be created, like `mkdir -p`; low risk).
8. **Help text.** Examples hard-code `~/tmp`, `-p 1000000` is shown as if typical, the first sentence of the `-d` note is unrelated to `-d`'s value, exit codes absent, `--norm`/`--plotend` explained with jargon ("1% read length intervals"), `-l` missing. README option list omits `--norm`.
9. **Messages.** Progress with `--verbose` only; a 100 GB BAM is silent for an hour, then "Writing to: ...". Add a one-line progress indicator to stderr when stderr is a TTY (reads processed, every ~1M reads), off otherwise. Warnings go through WARNING_MSG on stderr already (good). The "partial report" warning for `-p` is good, and the HTML should carry the same note (check; UNVERIFIED in this review).
10. **Subsampling semantics.** `-p N` = "first N reads". For sorted BAM that is a few kb of chr1. A random-ish sample needs either a full pass or index-based jumps; neither is needed for the common case (QC of everything). Keep `-p` with its honest name and document the bias; do not add `--sub`.

## 3. Proposed CLI for 2.4

Principle: `samstat file` always does the right thing. Everything else is optional and few.

```
samstat [options] <file> [<file> ...]       # or "-" for stdin
```

Zero required options. Final option set (8 flags, down from 11 + 2 dead):

| option | status | meaning |
|---|---|---|
| `-o, --outdir DIR` | rename of `-d/-dir` (old spelling kept) | output dir, created if missing; default: current directory. |
| `-p, --peek N` | keep | analyse first N reads. |
| `-t, --threads N` | keep, make real | decompression + parse threads. Default: min(4, ncores). |
| `--json` | keep | also write `<name>.samstat.json`. |
| `--multiqc` | keep | also write MultiQC `_mqc` files. |
| `--plotend` | keep for now | candidate for removal; see section 4. |
| `-q, --quiet` / `--verbose` | `--verbose` kept, `-q` added | `-q` silences the progress line. Optional; may skip. |
| `-h`, `-v/--version` | keep | |

Removed: `-l`, `--seed` (dead), `--norm` (becomes automatic behaviour for long reads; accepted and ignored).

Default output location: **change to current working directory?** Risk: breaks scripts that expect the report beside the input. Recommendation: keep "beside the input" as default (backwards compat, matches 2.3 docs and tests/goldens), but if that directory is not writable fall back to the CWD with a visible message, rather than failing at the end. For stdin the output is `stdin.samstat.html` in CWD unless `--name` ... (see below). Decision for the maintainer.

Stdin naming: `-` produces `samstat.html` (or `-o`-dir/`samstat.html`). No `--name` flag; users can rename. A second, optional `--prefix` is NOT proposed (complexity).

Input detection: content sniffing, extension as fallback for ambiguous cases (empty file). Accept `.cram` (needs reference for MD-based errors: use `REF_PATH`/`-T`; without reference, mismatch/MD profiles are unavailable. Decide whether CRAM is in scope; htslib 1.16 supports reading. UNVERIFIED how the current parser behaves without an MD tag - check in implementation).

Messages: keep the existing LOG/WARNING style, add (a) the detected read type ("long reads: median 8.2 kb, N50 14 kb - using long-read layout") to stderr and the HTML header, (b) a TTY progress line, (c) an exit-time summary line "wrote X.samstat.html (12.3M reads, 41 s)".

Exit codes: keep 0/1/2; print them in `-h`.

## 4. Remove / merge candidates (the "does it earn its place" list)

| item | verdict |
|---|---|
| `-l` | remove (dead). |
| `--seed` | remove (dead); reintroduce only with a real random subsample. |
| `--norm` | remove as flag; auto-on for long reads, otherwise absent. |
| `--plotend` | discuss; end-anchored plots matter for adapter/trimming QC in short reads. Either default-on for short reads or drop. Report-design question, not CLI. |
| `-p` | keep, it is cheap and useful for quick looks. |
| `-t` | keep but make it work for a single file. |
| `--verbose` | keep (debug). |
| `--json`/`--multiqc` | keep; consider folding into one `--format html,json,multiqc` only if a third machine format appears. Not now. |

## 5. Backwards compatibility

- Output filename `<input>.samstat.html`, `.samstat.json`, `.samstat_*_mqc.json` unchanged (golden tests, MultiQC globs, bioconda recipe).
- `-d`, `-dir`, `-peek` keep working (getopt_long_only already accepts them); add `-o/--outdir` as the documented name. Do not make `-d` an error.
- `--norm`, `-l`, `--seed` parse and are ignored with one deprecation line on stderr (exit code unchanged). Removal in 3.0.
- Exit codes unchanged.
- Test impact: tests with `-l`/`--seed` (if any; not enumerated in this review - grep tests/CMakeLists.txt) need updating; the `cli` ctest label covers option parsing and is where new cases go: stdin, magic sniffing without extension, read-only dir, name collision, created `-o` dir, `-t` on single file.

## 6. Priority order (smallest effort per researcher benefit)

1. Delete dead `-l`/`--seed` with deprecation (a few lines).
2. Up-front output-writable check and name-collision check; create `-o` dir.
3. Rewrite `-h` text: no hard-coded home paths, exit codes, one-line option lists, 3 examples (single file, many files with `-o`, pipeline with `--json --multiqc`). Update README option list.
4. Content sniffing + `-` for stdin (+ accept `.cram`).
5. Auto short/long layout with printed decision.
6. TTY progress line and final summary line.
7. Real `-t` for a single file (performance plan; CLI only needs the semantics).

## 7. Suggested `-h` (proposal)

```
samstat 2.4  QC for short and long reads (SAM/BAM/CRAM, FASTA/FASTQ)

Usage: samstat [options] <file> [<file> ...]      ("-" reads stdin)

  Writes <file>.samstat.html next to each input. No options are required.

Options:
  -o, --outdir DIR   write reports to DIR (created if needed)
  -p, --peek N       use only the first N reads            [all]
  -t, --threads N    threads                               [4]
      --json         also write <file>.samstat.json
      --multiqc      also write MultiQC files <file>.samstat_*_mqc.json
  -h, --help   -v, --version

Exit: 0 ok, 1 an input failed (others still processed), 2 usage error

Examples:
  samstat sample.bam
  samstat -o qc/ run/*.bam
  samtools view -h x.bam | samstat -
```
