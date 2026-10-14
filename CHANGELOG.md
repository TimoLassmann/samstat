# Changelog

All notable changes to SAMStat are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased] - 2.3.0

### Fixed

- Include paths for `plot.h` so the tree builds against the new tld layout.

### Added

- CTest suite (`ctest`) with smoke tests over every input in `dev/`.
- GitHub Actions: Ubuntu and macOS release builds, an AddressSanitizer job, and a tag-triggered release workflow publishing a source tarball (including the tld submodule) with its sha256.
- `CITATION.cff`, `CHANGELOG.md`, and a draft Bioconda recipe under `packaging/bioconda/`.

### Changed

- The `tld` library (`src/tld`) is now a submodule tracking upstream <https://github.com/timolassmann/tld> instead of a vendored copy.
- Minimum CMake version is 3.19.

### Removed

- Nothing yet.
