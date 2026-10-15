# Changelog

All notable changes to SAMStat are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [2.4.0] - unreleased

### Breaking changes

- **JSON format version 2.** `format_version` is 2; the version 1 keys are gone and the former `v2` block is the top level (`summary`, `mapq_bins`, `axes`, `read1`, plus `read2` for paired input). Per-base mismatch counts are replaced by one total per bin, and the quality-by-position matrix and clip-trimmed length no longer exist. Summary fractions are written with four decimals and are `null` when they do not apply. See `docs/json-format.md` and `docs/samstat.schema.json`. The MultiQC files change with it: seven general-stats columns, five sections, and a relative position axis for both position plots.
- **Protein FASTA is unsupported.** Inputs are now read with htslib, which treats a protein FASTA as nucleotides with every non-ACGT letter counted as N, so the report shows nearly 100% N. There is no detection and no error.
- **Unknown input exits 1 for that file.** Input type is detected from the content, not the extension. Non-empty input that is not SAM, BAM, FASTA or FASTQ fails that file with exit 1; the other files are still processed. Empty input is still a warning and exit 0.
- **`-p` counts.** With `-p N`, the secondary and supplementary counts cover only the reads up to the cutoff, so they no longer depend on batch size or `-t`.
- **Read length is the full read length** (soft-clipped bases included), not the clip-trimmed length. Lengths, N50 and base totals change for soft-clipped reads.
- **Deprecated options** `-l`, `--seed`, `--norm` and `--plotend` are accepted, print a warning and are ignored.

### Changed

- Intra-file threading: `-t` is now the total thread budget (default 4) and also speeds up a single input; output is byte-identical for any `-t`. With at least as many inputs as threads, inputs run one per thread as before.
- The HTML report is rebuilt: one page with a summary strip, optional notes, a MAPQ key table and six fixed plots as inline SVG; no Plotly, no JavaScript, no network request; prints cleanly to A4. The plots are read length, read mean quality, quality along the read, base composition along the read, error rates along the read (aligned input) and read identity (aligned input). Rates replace raw counts in the error profile.
- Reads are taken from SAM, BAM, FASTA and FASTQ through htslib; `-` reads standard input.
- Mismatches are taken from = and X operations in the CIGAR when a read has no MD tag; a read of = operations only counts as a perfect match. A note in the report states how many mapped reads lack MD. bwa-mem and bwa-mem2 write MD by default; STAR and minimap2 (without `--MD`) do not.

- `-t` accepts 1 to 64; larger values are rejected with exit 2.
- Memory is bounded by the thread count and the batch size, never by a setting or by the input: a record longer than 50 Mb (a genome FASTA, a corrupt record) fails that input with exit 1; a deletion in a CIGAR costs one alignment column however long it is; a buffer that held a long read is released after its batch; an unrecognised input is probed in fixed 64 KB chunks instead of a line at a time. The string buffer no longer copies its whole content on every append once it has grown (the `tld_strbuf_resize` fix). The test-only environment variables can only lower the built-in limits.
- Plot display: long-read lengths are drawn in 20 bins per decade (stored counts unchanged), the identity histogram of long reads uses the narrowest of 0.1/0.2/0.5/1 % bins that gives at most 50 bins, a MAPQ class is drawn in plots 3 and 5 only where at least 50 of its reads cover the position, plot 5 rate axes reach at least 0.5 %, count and position axes have whole-number ticks, and single figures use the full page width. Reads of one length get one sentence instead of a one-bar histogram (and a plot replaced by a sentence has no caption); the read-quality axis spans the data instead of starting at Q0; log length axes label 1, 2 and 5 of each decade when the labels fit.
- Every plot split by MAPQ class carries its own key above the panels, naming each class with its read count (and "too few to draw" when a class has reads but no line), so the reader no longer scrolls back to the MAPQ table.
- The per-position plots 3 to 5 are replaced by a one-line note for inputs with fewer than 50 reads.

- The tld git submodule is gone. The few parts of tld that samstat uses (error macros, allocation, string buffers, sequence buffers, file names, reverse complement) are in one header, `src/util/tldlite.h`, so a plain `git clone` and GitHub's own source archive build without extra steps. Output is byte-identical to the build with tld.

- Each input now prints one summary line on stderr (format, reads, N50, mean Q, % mapped, layout, time); a progress counter is shown only when stderr is a terminal.

### Added

- "Compare the MAPQ classes by shape": a collapsed second figure under plots 1, 2 and 6 of aligned input, with each MAPQ class scaled to 100 % of its own reads, so that a small class (unmapped, low MAPQ) can be compared with the bulk by shape. Classes with fewer than 50 reads are left out; the figure is shown only when at least two classes remain.
- CRAM input. The reference is found through `REF_PATH` or `REF_CACHE` only (no network access); a CRAM whose reference cannot be found fails that file with exit 1 and a message naming both variables. The bundled htslib is built with bzip2 and xz support.

### Removed

- Plotly and all JavaScript from the HTML report.

## [2.3.0] - 2026-10-14

### Fixed

- **Soft-clipped reads:** a soft clip at the 5' end (and a soft clip next to a hard clip) was handled wrongly, so base composition, quality and mismatch positions were shifted and read lengths were too long (for example `5H10S40M` was counted as 50 bases, now 40). Numbers for soft-clipped reads change.
- **Supplementary alignments** (flag 0x800) are no longer counted as reads, like secondary ones. Read totals and all per-read statistics for such files change. Both kinds are now counted separately (`secondary` and `supplementary` in the JSON summary and the HTML overview).
- FASTQ base qualities were decoded with Phred offset 0 instead of 33, which put every quality 33 too high in the quality statistics (and the JSON and MultiQC files); they are now Phred+33.
- FASTA/FASTQ inputs without reads (an empty `.fq.gz`, only newlines, header-only FASTA) write an empty report with a warning and exit 0; a last record without trailing newline is no longer dropped or rejected.
- `meta.samstat_version` in the JSON report was the htslib version; it is SAMStat's.
- MultiQC sample names strip the known extensions case-insensitively (`READS.FQ` gives `READS`); `.cram` is no longer in the list because CRAM is not supported.
- Reads with missing quality (`*`) no longer count as Q41 in `--plotend` and `--norm` plots; base qualities are clamped to 0-41.
- Crashes and memory errors: heap overflow for reads of 1024*n bases, records with SEQ `*`, out-of-bounds MD tag parsing, bad indexing for histograms and plots, missing format arguments in log messages.
- Unreadable, non-BAM or truncated inputs now give a clear error instead of a crash or an empty report.
- Failed inputs no longer stop the run: the other inputs are still processed, and samstat exits 1 at the end. Usage errors exit 2; `-h` and `-v` exit 0.
- File type detection ignores case and recognises more FASTA/FASTQ suffixes; an unknown extension is a usage error.
- `-t`, `-p` and `-l` are validated strictly (digits only, overflow checked); the `-l` upper bound check was unreachable.
- `-p` now analyses exactly N reads (previously rounded up to a whole chunk); the minimum is 1.
- Empty inputs (no reads) write a report with a warning and exit 0; reports no longer contain `nan` or `inf`.
- Output path length and write errors are checked when writing reports.
- HTML report: valid CSS, `lang`, charset and viewport, escaped file names. Plotly is pinned to 1.58.5 instead of `latest` (the unused d3 is no longer loaded) and line traces use `scatter` with `mode: lines`.
- The footer cites SAMStat 2 (Lassmann, Bioinformatics 2023).
- Resources are freed on all error paths.
- Include paths for `plot.h` so the tree builds against the new tld layout.

### Added

- `--multiqc` option: writes MultiQC custom content files `<input>.samstat_{general,mapq,length,basequal}_mqc.json` (see `docs/multiqc.md`).
- `--json` option: writes a machine-readable `<input>.samstat.json` next to the HTML report (format in `docs/json-format.md`, schema in `docs/samstat.schema.json`).
- `SAMSTAT_SANITIZE` CMake option (AddressSanitizer and UBSan) and `SAMSTAT_USE_SYSTEM_HTSLIB` option.
- CTest areas: `build`, `json`, `packaging`, and `golden` (reports of the 2.3.0 output; differences from the 2.2 baseline are explained in `docs/golden-diff-v2.3.md`; regenerate with the `update-goldens` target). Hand-written edge-case fixtures in `tests/data/`.
- CTest suite (`ctest`) with smoke tests over every input in `dev/`.
- GitHub Actions: Ubuntu and macOS release builds, an AddressSanitizer job, and a tag-triggered release workflow publishing a source tarball (including the tld submodule) with its sha256.
- `CITATION.cff`, `CHANGELOG.md`, and a draft Bioconda recipe under `packaging/bioconda/`.

### Changed

- The `tld` library (`src/tld`) is now a submodule tracking upstream <https://github.com/timolassmann/tld> instead of a vendored copy.
- HDF5 is no longer a dependency: tld (961795c) makes its hdf5 wrapper optional and samstat builds it with `TLD_USE_HDF5=OFF`.
- Minimum CMake version is 3.19.
- `config.h` is generated in the build tree; the configure step no longer updates submodules; htslib is built with proper CFLAGS and without autoreconf.

### Removed

- Dead metrics and report code (`src/metrics/`, `src/report/report.c`) that was not part of the build.
