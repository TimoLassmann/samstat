# Topic: ecosystem (FINAL design for 2.4)

Base: the **minimal** proposal. Grafted from the researcher proposal: one MultiQC error-profile linegraph (so SAMStat's distinctive plot appears in the cohort view), log-spaced length bins in the MultiQC length file, a schema-validation test over all goldens, and a local bioconda-utils dry run before the PR.

Principle: users get one HTML report per sample and use MultiQC for a cohort. Every ecosystem item is either a way to distribute SAMStat or a view of the one versioned JSON contract. No integration layer computes new metrics.

## 1. Decisions

| Item | Decision |
|---|---|
| Release tarball, bioconda, biocontainers | **IN, first.** It is the gate for everything else. |
| JSON `format_version` 2 as the single contract | **IN.** One breaking bump in 2.4.0, then frozen for 2.x. |
| `--multiqc` custom content | **KEEP.** Add two general-stats columns and one error linegraph, and switch the length file to log bins. |
| nf-core module | **IN, after the bioconda 2.4.0 package exists.** It lives in nf-core/modules, not in this repo. |
| Native MultiQC module | **OUT for 2.4.** Open one module-request issue after bioconda. |
| Galaxy wrapper | **OUT.** Accept a community PR if one appears. |
| bio.tools entry | **IN, 10 minutes.** Whether an entry already exists is UNVERIFIED. |

## 2. What already exists (verified in code, 2026-10-14)

- `packaging/bioconda/meta.yaml` is a draft. It has version 2.3.0, a placeholder sha256 of zeros, and a test that runs `samstat -v` only.
- `packaging/bioconda/build.sh` **already exists**. It runs cmake with `-DGIT_SUBMODULE=OFF` and the bundled htslib. Neither proposal noticed it.
- `packaging/make_source_tarball.sh` exists. `.github/workflows/release.yml` checks out with `submodules: recursive` and attaches `dist/samstat-*.tar.gz` and its `.sha256`. **The tarball work is done.** Only the recipe needs updating.
- `tests/cmake/check_packaging.cmake` (`packaging:metadata`) already asserts that the recipe version equals `PROJECT_VERSION`, and checks CHANGELOG and CITATION.cff.
- `samstat -v` prints `SAMStat 2.3.0`. The code is `param.c:125` `"%s %s\n"` with `PACKAGE_NAME` = `CMAKE_PROJECT_NAME`, and `project(SAMStat ...)` is set in `CMakeLists.txt:3`. The researcher proposal's lower-case `samstat 2.3.0` is wrong. The CLI tests (`tests/cli/CMakeLists.txt:80-81`) only check `[0-9]+\.[0-9]+`.
- The JSON has `"format_version": 1` at `json_report.c:75`. `docs/json-format.md:13` already defines the policy: a breaking change bumps the version, a new optional field does not.
- The MultiQC general-stats column ids, from `mqc_report.c:179-193`, are `samstat_total_reads`, `samstat_pct_mapq_ge30`, `samstat_pct_unmapped`, and `samstat_pct_proper_pair` (paired data only). The sections are general, mapq, length and basequal. There is no error profile. File names follow `<input>.samstat_<section>_mqc.json`.
- `-d` must name an existing directory (`param.c:199-208`). Adding `-o/--outdir` belongs to the CLI topic.

## 3. The contract (must be fixed before collector and report work starts)

1. **JSON `format_version` 2.**
   - The 2.4 layout changes are breaking: bounded or relative position bins, a log length histogram with N50, an identity histogram, and a raised quality cap.
   - Bump to 2 once, in 2.4.0.
   - Update `docs/json-format.md` and `docs/samstat.schema.json` in the same commit.
   - Do not emit version 1 in parallel. The changelog gets a one-line note for version 1 consumers.
   - After 2.4.0, only optional fields may be added in the 2.x line.
   - **The collector and report topics own the field layout.** This topic requires three things of it:
     - the field names are reserved before implementation starts;
     - the JSON size is independent of read length;
     - every existing stable id (MAPQ bin ids `mapq_ge30`, `mapq_10_29`, `mapq_lt10`, `unmapped`) keeps its name.
2. **`samstat -v` stays `SAMStat X.Y.Z`** on stdout with exit 0.
   - The nf-core module strips the prefix with `sed 's/^SAMStat //'`.
   - Changing the output would only benefit that one consumer, so it stays as it is.
   - Tighten the CLI test to the full regex `^SAMStat [0-9]+\.[0-9]+\.[0-9]+$`.
3. **Output names are deterministic from the input name.** They are `<input>.html`, `<input>.json` and `<input>.samstat_<section>_mqc.json`. No timestamps and no absolute paths appear in any output. This is already true per `docs/multiqc.md`, and it must stay true.
4. **The MultiQC column ids are a contract.**
   - The existing four ids are never renamed.
   - New ids: `samstat_n50`, plus `samstat_median_identity` for aligned input only. Neither column is written when it is not applicable.

## 4. MultiQC custom content: the 2.4 changes (about 1 day)

The aim is that a cohort reader can answer "which sample is the outlier, and is it the reads or the mapping?" The output must stay small enough for 5000 samples. MultiQC's scaling requirement is at https://github.com/MultiQC/MultiQC/blob/main/.github/CONTRIBUTING.md.

- **General stats.** Keep the existing four columns. Add `samstat_n50` and `samstat_median_identity`.
- **Length file.** Use the same log-spaced bins as the HTML and JSON. This replaces the current linear series, which is block-summed above 2000 points.
- **New `samstat_error` linegraph.**
  - One line per sample, showing the combined mismatch + insertion + deletion rate per aligned base, for MAPQ >= 30 reads only.
  - It uses the same bounded position axis as the report (relative bins for long reads), so it has a fixed small number of points.
  - It is written only for aligned input.
  - This is the single place SAMStat's differentiator reaches the cohort view.
- **Not in MultiQC:** per-MAPQ-bin plots, base composition, read 2 and identity histograms. All of these stay in the per-sample HTML.
- **Sample-name cleaning.** Keep the documented `fn_clean_sample_names: false` workaround (`docs/multiqc.md`). Move it to the top of the page. No code change.

## 5. Bioconda (0.5 day plus review wait)

- In `packaging/bioconda/meta.yaml`:
  - set the version to 2.4.0 and the real sha256 from the release asset;
  - add `extra: identifiers: [biotools:samstat]`;
  - extend `test: commands` with one real run, `samstat --json --multiqc` on a tiny SAM. Write the SAM inline with `printf` in the test, because the recipe's ability to ship test source files is UNVERIFIED.
- `build.sh` stays as it is. The bundled htslib is built offline from `thirdparty/`.
- **Before opening the PR,** build the recipe locally with `bioconda-utils build` (https://bioconda.github.io/contributor/workflow.html, UNVERIFIED in detail). The risk this addresses is that the ExternalProject htslib configure step needs autoconf, libcurl or network access inside the bioconda sandbox. Fall back to the conda `htslib` package only if this build fails.
- Never put release candidates in bioconda.
- Biocontainers images are built automatically from merged bioconda recipes, so there is no separate work. This is UNVERIFIED in detail.
- This also resolves issues #6 (Conda availability) and #13 (build failures, in part).

## 6. nf-core module (about 1 day, after bioconda 2.4.0)

- Use `modules/nf-core/cramino` as the template: https://github.com/nf-core/modules/tree/master/modules/nf-core/cramino. Guidelines are at https://nf-co.re/docs/guidelines/components/modules.
- The module consists of `main.nf`, `meta.yml`, `environment.yml` (pinned `bioconda::samstat=2.4.0`) and `tests/main.nf.test` with a `-stub` test.
  - Outputs: `*.html`, `*.json` and `*_mqc.json`.
  - Versions go through `topic: versions`, not `versions.yml` (per `research/multiqc-nfcore.md`).
  - Snapshots are stable because the outputs are deterministic.
- Run it with `-d .` on the input staged as `${prefix}.bam`.
- Label it `process_single`. Once single-file threading ships (perf topic), switch it to `process_low` with `-t ${task.cpus}`.
- Test data comes from existing nf-core/test-datasets. Which files to use is UNVERIFIED.
- No copy of the module is kept in this repo.
- The pitch, given that fastqc, nanoplot and cramino are already in the module list, is the MAPQ-stratified error profile, which none of those tools offer.

## 7. Offline HTML: a release gate, owned elsewhere

The presentation topic owns the fix (inline SVG). The ecosystem requirement is that **the report makes zero network requests**, because nf-core runs on clusters without internet.

## 8. Explicitly OUT

- A native MultiQC module in 2.4. It would need 2-3 days of Python, CI (`--strict`, mypy, ruff), a test-data PR and an unbounded review wait (https://docs.seqera.io/multiqc/development/modules). Triage scores conda downloads, which are zero before bioconda (https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md).
- Deleting `mqc_report.c`. Do this only if a native module ships, and then only after one deprecation cycle.
- A Galaxy wrapper, which needs a version-bump PR on every release (https://github.com/galaxyproject/tools-iuc).
- Switching to the conda `htslib` package, unless the local bioconda build fails. The htslib version bump is owned by the perf topic.
- Supporting JSON versions 1 and 2 in parallel, a second JSON profile, or any TSV or Parquet output.
- A MultiQC config file shipped with SAMStat.
- A CI job that installs into a clean conda environment. The local `bioconda-utils` build covers the same risk once per release.

## 9. Order and effort

| Step | Effort | Depends on |
|---|---|---|
| Reserve `format_version` 2 field names and the MultiQC column ids | 0.5 d | before collector and report implementation |
| Implement the 2.4 features (other topics) | - | - |
| `--multiqc` changes (§4) | 1 d | JSON layout fixed |
| Recipe update and local `bioconda-utils` build | 0.5 d | code complete, htslib bump done |
| Tag v2.4.0; release workflow attaches tarball and sha256 (already automated) | 0.1 d | above |
| Bioconda PR | 0.2 d plus wait | tag |
| bio.tools entry | 0.1 d | none |
| nf-core module | 1 d | bioconda 2.4.0 merged |
| MultiQC module-request issue with example JSON | 0.2 d | bioconda 2.4.0 merged |
| **Total** | **about 3.6 d** | |

## 10. Risks

- **The bioconda build fails on the bundled htslib** (ExternalProject, autoconf, libcurl). Mitigation: the local `bioconda-utils` build. Fallback: the conda `htslib` package.
- **The JSON layout changes after the nf-core snapshot is taken.** Mitigation: the freeze rule in §3 and the schema test.
- **Review latency** at bioconda and nf-core is outside our control. Neither one blocks the GitHub release.
- **MultiQC behaviour changes.** `docs/multiqc.md` records the validated version (1.35). Re-run `multiqc:run` at each release.

## 11. Tests

- `cli`:
  - `samstat -v` matches `^SAMStat [0-9]+\.[0-9]+\.[0-9]+$` and exits 0. This tightens `version_short` and `version_long`.
  - With `-d <dir>` and a symlinked input, the output names are `<input-name>.*`, and no output contains the absolute input path.
- `packaging`: the existing `packaging:metadata` already checks the recipe version. Add a check that the recipe test block runs more than `samstat -v`. The sha256 placeholder check runs at release time only, in `release.yml`, not in CTest.
- `json`:
  - `format_version == 2`.
  - Every JSON golden validates against `docs/samstat.schema.json`.
  - A size bound for a 1 Mb read fixture. That test is owned by the collector topic.
- `multiqc`:
  - The exact general-stats id list is pinned in a golden.
  - Hand-computed values for `samstat_n50`, `samstat_median_identity` and the `samstat_error` linegraph.
  - The `samstat_error` file is absent for FASTQ input.
  - `multiqc:run` asserts that the new section appears for all three samples.
- `report`: the HTML contains no `http://` or `https://` script, stylesheet or font references.
- External: `nf-core modules lint samstat` and `nf-core modules test samstat`, including the stub.

## 12. DECISION NEEDED

1. **JSON break.** Confirm a clean `format_version` 2 with no version 1 output, rather than keeping version 1 and adding the new data as optional fields only. The binned position axis cannot be expressed additively.
2. **MultiQC error linegraph.** Include it in 2.4 (+0.5 d, recommended), or ship only the new general-stats columns.
3. **Native MultiQC module.** Are you willing to maintain a Python module upstream at all? If not, the module-request issue should say so, and custom content becomes permanent.
4. **nf-core module authorship.** Do you submit and maintain it under your GitHub account, or hand it to an nf-core volunteer? Each version bump needs someone.
5. **Release timing.** Should bioconda wait for the htslib 1.24 bump (perf topic), or ship 2.4.0 on 1.16 if that bump slips?
