# Topic: release-cleanup (FINAL)

The minimal proposal is the base. From the researcher proposal I took four things: the offline clean-build gate before tagging, the edge-case FASTQ fixtures that must exist before the reader is swapped, the rule that htslib is bumped before the performance work, and the CRAM no-network test. I checked the claims about the code against the repo on branch `planning/v2.4` (HEAD c139e9d), and five corrections to both proposals follow.

## 0. Repo facts this design relies on (verified)

| Claim | Status in repo |
|---|---|
| 2.3.0 is untagged | Only `v2.0.0` exists. `CHANGELOG.md` has `## [2.3.0] - unreleased`. |
| Release tarball contains the tld submodule | **Already true.** `.github/workflows/release.yml` checks out with `submodules: recursive` and runs `packaging/make_source_tarball.sh`, which archives each submodule at its pinned commit. |
| Release tarball contains htslib | **Already true.** `thirdparty/htslib-1.16.tar.bz2` is git-tracked, so `git archive` includes it. |
| Bioconda recipe | `packaging/bioconda/meta.yaml` has a placeholder sha256. **`build.sh` already exists** (cmake configure, build, install, `-DGIT_SUBMODULE=OFF`). The only test is `samstat -v`. |
| Recipe host deps | Lists `zlib bzip2 xz libcurl`, but `thirdparty/htslib.cmake` configures htslib with `--disable-bz2 --disable-lzma --disable-libcurl --disable-s3 --without-libdeflate`. Only zlib is used. |
| htslib pin | `htslib_URL_MD5` (MD5, not SHA256). The system-htslib option needs `>=1.16` (`CMakeLists.txt:51`, `thirdparty/htslib.cmake:4`). |
| `samstat -v` | Prints `SAMStat 2.3.0` (`param.c:125`), **not** a bare version. |
| CHANGELOG | **`CHANGELOG.md` already exists** (Keep a Changelog format). The legacy `ChangeLog` and `NEWS` files also exist and stop at 2.2.3. |
| tld strbuf bug | `tld_strbuf_resize` (`src/tld/src/string/*.c:40`) swaps the buffer but never sets `b->alloc_len = new_size`, so every append past the first allocation reallocates. The fix is one line. |
| JSON workaround | `json_report.c:31` sets `JSON_BUF_BYTES 65536`, with `json_flush` calls at lines 127, 147 and 158. |
| FASTA/FASTQ workaround | `prepare_fasta_fastq_input` (`samstat.c:257`, called at line 344). |
| Plotly | Loaded from the CDN at 1.58.5 (`stat_report.c`). |
| Dead options | `-l` sets `p->report_max_len`, which is read only inside a comment (`collect.c:594`). `--seed` sets `p->seed`, which is read nowhere. |
| Version source | `project(SAMStat VERSION "2.3.0")` in `CMakeLists.txt:3`. The `packaging:metadata` and `packaging:tld_pin` tests exist. |

## 1. Decisions, in execution order

### R1. Tag v2.3.0 as it is (0.25 day)
- Push `main` and get every CI job green (ubuntu, macOS, ASAN, multiqc). Change no code except CI fixes.
- Change `## [2.3.0] - unreleased` to the release date in `CHANGELOG.md` (one commit).
- **Gate before tagging:** run `packaging/make_source_tarball.sh 2.3.0 dist`, extract the tarball in a clean temp directory, and run the build with networking disabled. Then run `samstat -v` and `ctest -L "build|cli"`. The workflow probably already works, but this has never been exercised offline.
- Tag `v2.3.0` on main. The existing release workflow attaches the tarball and its `.sha256`.
- Do not hotfix 2.3.x for long-read memory, the single core per file or the Q41 cap. Those are 2.4 work.
- Keep `samstat -v` as `SAMStat 2.3.0`. nf-core and bioconda parse it with `sed`/regex, and changing it gains nothing.

### R2. Submit bioconda for 2.3.0 now (0.5 day plus review wait)
- In `meta.yaml`:
  - Use the real sha256 from the release asset.
  - Drop the `bzip2`, `xz` and `libcurl` host deps, which are unused (see §0).
  - Remove the "DRAFT" header.
  - Keep the bundled htslib. Switch to conda `htslib` with `-DSAMSTAT_USE_SYSTEM_HTSLIB=ON` only if bioconda CI rejects the ExternalProject build.
- Add a second test command that runs `samstat` on a 3-read FASTQ written inline with `printf` in `test.commands`. This avoids depending on `test.source_files` (UNVERIFIED whether that is allowed for URL sources).
- Open the PR to https://github.com/bioconda/bioconda-recipes. No samstat recipe exists today, per research/multiqc-nfcore.md; I did not re-check this (UNVERIFIED).
- Once the recipe is merged, the bioconda autobump bot handles 2.4.0. This is the reason to do it now rather than wait for 2.4.

### R3. Fix the tld strbuf bug upstream and bump the pin once (0.25 day)
- Fix `tld_strbuf_resize` in https://github.com/timolassmann/tld by adding `b->alloc_len = new_size;`. Add an upstream test that does 1e6 appends and asserts `alloc_len >= len` and fewer than 64 reallocations.
- In samstat:
  - Bump the `src/tld` pin.
  - Update `tests/cmake/check_tld_pin.cmake`.
  - Delete `JSON_BUF_BYTES` and the periodic `json_flush` calls in `json_report.c`. Build the whole JSON in one strbuf and write it once. The 2.4 bounded-bin collector keeps the JSON small, so streaming is no longer needed.
- Bump the tld pin **once** in 2.4, together with any other tld change the 2.4 topics need. Do not bump it per fix.

### R4. The tld FASTA/FASTQ reader: retire it rather than fix it
- The I/O topic replaces the tld reader with htslib `sam_open`/`sam_read1`. Per research/code-io.md, htslib 1.16 already handles unterminated last records, `.gz` and pipes, and reads 1 M reads in 0.78 s gzipped against 1.78 s today.
- That deletes `prepare_fasta_fastq_input` (the double gunzip and the temp copy) and enables stdin. The I/O topic owns this work. Release-cleanup only requires the protective fixtures listed under §4 to land **before** the swap.
- Fix the unterminated-record bug in tld upstream only if the tld reader stays (see DECISION NEEDED 1).
- Cost: protein FASTA support goes away, because htslib stores nucleotides only. The change needs a CHANGELOG line.
- Open check for the I/O topic: htslib's FASTQ R1/R2 handling (`/1` `/2` suffixes, CASAVA flags) against the current `fastx` and `golden` tests. UNVERIFIED.

### R5. Bump bundled htslib from 1.16 to the current release (target 1.24) in an isolated commit, before the performance work (0.5 to 1 day)
- **Why:** per https://github.com/samtools/htslib/blob/develop/NEWS, as summarised in research/perf-htslib.md and not re-fetched here (UNVERIFIED):
  - CRAM CVE fixes in 1.20 and 1.23.1.
  - Faster `bgzf_read` in 1.21.
  - No default EBI reference fetch, plus pipe and threaded-EOF fixes, in 1.22.
  - MN tag support in 1.18, relevant to the modbase topic.
  - 1.24 is listed as released 9 Jul 2026 (UNVERIFIED).
- **Mechanics:**
  - Replace `thirdparty/htslib-1.16.tar.bz2`.
  - Switch `URL_MD5` to `URL_HASH SHA256=`.
  - Keep the configure flags unchanged in this commit.
- **Order:** do this commit first, alone. Run the full suite and the goldens. Any golden diff is explained in the commit message before anything is regenerated. This gives the performance work a final baseline.
- **Constraint for the performance topic:** from 1.23.1, setting a thread pool twice on a file is disallowed. Create one shared `hts_tpool` per run and attach it once per opened file.
- **Constraint for the I/O topic:** if CRAM input is accepted, the I/O commit (not this one) re-enables `bz2` and `lzma` in the configure line, because CRAM blocks may use those codecs. bzip2 and xz then go back into the recipe host deps.
- **Not in scope:** keep `--without-libdeflate`, `--disable-libcurl` and `--disable-s3`. libdeflate is a possible speed-up for the performance topic to benchmark (UNVERIFIED gain). It is not a release item.
- The `SAMSTAT_USE_SYSTEM_HTSLIB` minimum stays `>=1.16` unless a 2.4 feature needs a newer API (see DECISION NEEDED 2).

### R6. Plotly: replace it, do not upgrade it
- Do not upgrade Plotly. Inlining it costs 1.2 MB for basic, 1.51 MB for cartesian and 4.83 MB for full, per https://github.com/plotly/plotly.js/blob/master/dist/README.md (verified in research/presentation.md). The CDN breaks offline, HPC and container use.
- The report topic owns the replacement. Recommended: inline SVG generated in C. Fallback: an inlined, pinned plotly-basic. Never a CDN.
- Release-cleanup adds only the guard test in §4. No interim vendored Plotly ships in 2.3.x.

### R7. Dead options (0.1 day; shared with the CLI topic)
- In 2.4, `-l/--len` and `--seed` still parse but print one stderr line, `samstat: warning: -l is ignored and will be removed in 2.5`, and exit 0.
- Remove them from `print_help` and remove `report_max_len` and `seed` from `param.h`.
- Remove both options entirely in 2.5.

### R8. One changelog, one version source (0.1 day)
- `CHANGELOG.md` is the only changelog. Add a single line at the top of the legacy `ChangeLog` and `NEWS`: "Superseded by CHANGELOG.md". Do not delete them, because old tarballs reference them.
- The version lives only in `project(VERSION)`. The release workflow already derives the tarball name from the tag. Add one check to `packaging:metadata` that the tag-less `project(VERSION)` matches the top CHANGELOG heading that is not "unreleased" at release time. Skip this check if it turns out to be fiddly.

### R9. Release 2.4.0
Use the same workflow and the same offline-tarball gate as R1. The bioconda bot opens the version bump, and it needs only a review.

## 2. Explicitly OUT
- Galaxy tool, nf-core module, native MultiQC module. Keep the 2.3 `_mqc` output, and revisit after bioconda has downloads.
- A Docker or Singularity image (BioContainers builds one from bioconda automatically), Homebrew, Windows.
- A Plotly upgrade or a vendored Plotly in 2.3.x.
- Making libdeflate, libcurl or S3 required.
- Changing the `samstat -v` format.
- A tld rewrite. The only tld change is the strbuf fix, plus the reader fix only if the reader stays.
- Release automation beyond the existing tag workflow.

## 3. Affected files

| File | Change | Owner |
|---|---|---|
| `CHANGELOG.md`, `ChangeLog`, `NEWS` | date 2.3.0; supersede notes | release |
| `packaging/bioconda/meta.yaml` | sha256, host deps, extra test, drop DRAFT | release |
| `src/tld` (submodule), `tests/cmake/check_tld_pin.cmake` | pin bump | release |
| `src/report/json_report.c` | remove 64 KB flush workaround | release |
| `thirdparty/htslib-*.tar.bz2`, `thirdparty/htslib.cmake` | 1.24 tarball, SHA256 hash | release |
| `src/param/param.c`, `src/param/param.h` | `-l`/`--seed` deprecation | release/CLI |
| `tests/packaging/`, `tests/report/` | tests in §4 | release (via tests owner) |
| `src/samstat.c` (`prepare_fasta_fastq_input`), `src/htsinterface/htsglue.c` | reader swap | I/O topic |
| `src/report/stat_report.c` | Plotly removal | report topic |

## 4. Tests

| Test | Label | Lands with |
|---|---|---|
| Offline tarball build: run `make_source_tarball.sh`, extract to a temp dir, configure and build with no git and no network, check that `samstat -v` matches `^SAMStat [0-9]+\.[0-9]+\.[0-9]+$`. Run in CI as part of the release workflow (and optionally a nightly job), not in plain ctest, because it needs git. | packaging | R1 |
| The existing suite and goldens are unchanged after the htslib bump alone. | all | R5 |
| A CRAM file whose reference is unavailable fails fast with a clear message and exit 1, and makes no network access (`REF_PATH`/`REF_CACHE` unset, run under a timeout). | reader | R5, if CRAM is accepted by the I/O topic |
| FASTQ and FASTA with no trailing newline, each plain and `.gz`; an empty `.gz`; multi-line FASTA. Expected read counts come from the fixture definition, not from samstat output. Most already exist in `tests/fastx`; add only what is missing. | fastx | before R4 |
| `-l 5` and `--seed 3` exit 0 and print exactly one warning line on stderr. | cli | R7 |
| The generated HTML contains no `<script src="http` or `<link href="http`, and no `cdn.plot.ly`. | report | when the report topic lands |
| The JSON report for an existing fixture is byte-identical before and after the strbuf/flush change. | json | R3 |
| `check_tld_pin` updated to the new pin. | packaging | R3 |

## 5. Risks
- **Bioconda rejects the ExternalProject build** (network or autotools in the sandbox). Fallback: `-DSAMSTAT_USE_SYSTEM_HTSLIB=ON` with conda `htslib`. This works only if the system-htslib code path is CI-tested. Add one CI matrix entry with `SAMSTAT_USE_SYSTEM_HTSLIB=ON` (apt `libhts-dev`) before submitting.
- **htslib bump shifts output.** It is isolated in its own commit and diffed against the goldens before any other 2.4 change, so a regression has one source.
- **Dropping protein FASTA** is a user-visible regression. It needs a CHANGELOG entry and maintainer sign-off.
- **The bz2/lzma codecs are disabled** in the bundled htslib. CRAM files that use them will fail to decode until the I/O topic re-enables them.
- **UNVERIFIED:** the htslib 1.24 release date and NEWS contents (taken from the research notes), that bioconda has no recipe, the libdeflate gain, and htslib FASTQ R1/R2 semantics.

## 6. Effort
About 2 days of work plus the bioconda review wait. Critical path: R1, then R2, then R5 (before the performance work), then R3. R4 and R6 land inside the I/O and report topics.

## DECISION NEEDED
1. **Drop protein FASTA** (retire the tld FASTA/FASTQ reader for htslib)? Recommended: yes. If not, fix the unterminated-record bug in tld upstream and keep `prepare_fasta_fastq_input`.
2. **htslib target and system minimum.** Bundle 1.24, or the latest release at implementation time? Should the `SAMSTAT_USE_SYSTEM_HTSLIB` minimum be raised (for example to 1.18 for MN, if the modbase block is accepted)? Recommended: bundle the latest. Keep the minimum at 1.16 unless a feature needs more.
3. **Plotly replacement** (shared with the report topic): inline C-generated SVG (recommended), or an inlined plotly-basic of about 1.2 MB.
4. **Bioconda timing:** submit 2.3.0 now (recommended), or wait for 2.4.0.
