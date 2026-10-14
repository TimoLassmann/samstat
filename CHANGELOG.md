# Changelog

All notable changes to SAMStat are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
