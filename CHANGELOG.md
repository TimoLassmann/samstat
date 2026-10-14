# Changelog

All notable changes to SAMStat are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased] - 2.3.0

### Fixed

- Include paths for `plot.h` so the tree builds against the new tld layout.

### Added

- `--json` option: writes a machine-readable `<input>.samstat.json` next to the HTML report (format in `docs/json-format.md`, schema in `docs/samstat.schema.json`).
- `SAMSTAT_SANITIZE` CMake option (AddressSanitizer and UBSan) and `SAMSTAT_USE_SYSTEM_HTSLIB` option.
- CTest areas: `build`, `json`, `packaging`, and `golden` (baseline reports of the pre-fix behaviour; regenerate with the `update-goldens` target). Hand-written edge-case fixtures in `tests/data/`.
- `docs/tld-hdf5-followup.md`: upstream tld patch to make HDF5 optional.
- CTest suite (`ctest`) with smoke tests over every input in `dev/`.
- GitHub Actions: Ubuntu and macOS release builds, an AddressSanitizer job, and a tag-triggered release workflow publishing a source tarball (including the tld submodule) with its sha256.
- `CITATION.cff`, `CHANGELOG.md`, and a draft Bioconda recipe under `packaging/bioconda/`.

### Changed

- The `tld` library (`src/tld`) is now a submodule tracking upstream <https://github.com/timolassmann/tld> instead of a vendored copy.
- Minimum CMake version is 3.19.
- `config.h` is generated in the build tree; the configure step no longer updates submodules; htslib is built with proper CFLAGS and without autoreconf.

### Removed

- Dead metrics and report code (`src/metrics/`, `src/report/report.c`) that was not part of the build.
