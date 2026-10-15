# Releasing SAMStat and getting it into Bioconda, nf-core and MultiQC

Everything here is prepared and tested locally. The steps that publish something are yours.

| What | Where | What it gives users |
|---|---|---|
| Bioconda | `packaging/bioconda/` | `conda install -c bioconda samstat`, plus Docker and Singularity images (BioContainers) built automatically |
| nf-core module | `packaging/nf-core/samstat/` | any Nextflow pipeline can add SAMStat with `nf-core modules install samstat` |
| MultiQC module | `packaging/multiqc/` | SAMStat appears as a proper section in every MultiQC report, read from `samstat --json` |

Order: release, then Bioconda, then nf-core (it needs the Bioconda container). MultiQC is independent
and can go at any time.

## 1. Release (about 30 minutes)

1. Push both branches: `git push origin develop/v2.4 main`. Wait until the GitHub Actions checks are green
   (Ubuntu, macOS, sanitizer and MultiQC jobs).
2. Tag and push: `git tag -a v2.4.0 -m "SAMStat 2.4.0" && git push origin v2.4.0`. The release workflow
   runs the tests and creates the GitHub release with generated notes.
3. Checksum of the source archive (GitHub makes it; no submodules, so it is complete):

   ```sh
   curl -sL https://github.com/TimoLassmann/samstat/archive/refs/tags/v2.4.0.tar.gz | shasum -a 256
   ```

## 2. Bioconda (about 1 hour, then review)

1. Fork https://github.com/bioconda/bioconda-recipes and create a branch.
2. Copy `packaging/bioconda/meta.yaml` and `build.sh` to `recipes/samstat/`. Put the checksum from step
   1.3 in place of the zeros on the `sha256:` line.
3. Open a pull request. Bioconda's CI builds Linux and macOS (x86 and ARM) and runs the test commands of
   the recipe (they were run locally on macOS: all pass). A volunteer reviews; when it is merged the
   package and the container `quay.io/biocontainers/samstat:2.4.0--<build>` follow automatically.
4. Later versions: a bot opens the update pull request when you tag a new release.

## 3. nf-core module (about 1 hour, then review)

Needs the Bioconda package and its container from step 2.

1. Install the tools once: `pip install nf-core` (or `uv tool install nf-core`), Nextflow, nf-test,
   and Docker or Conda.
2. Fork https://github.com/nf-core/modules, create a branch, and copy `packaging/nf-core/samstat/` to
   `modules/nf-core/samstat/`.
3. In `main.nf`, replace the two `FIXME` container tags with the tag shown on
   https://quay.io/repository/biocontainers/samstat?tab=tags.
4. Generate the test snapshot and lint:

   ```sh
   nf-core modules test samstat --profile docker --update
   nf-core modules lint samstat
   ```

5. Open a pull request. nf-core reviewers usually ask for a round of small changes. Questions:
   #modules channel on the nf-core Slack.

The module runs `samstat --json -t <cpus>` and emits `html`, `json` (for MultiQC) and, with
`--multiqc` in `ext.args`, the `mqc` custom-content files. The three test inputs from nf-core's
test data (paired Illumina BAM, PacBio BAM, Illumina FASTQ.gz) were run locally with that command.

## 4. MultiQC module (about 1 hour, then review)

`packaging/multiqc/0001-New-module-SAMStat.patch` adds the module to MultiQC
(`multiqc/modules/samstat/`, registration, 7 unit tests). It reads `*.samstat.json` (format version 2)
and shows: General Statistics (reads, MAPQ >= 30, unmapped, N50, mean Q, identity), reads per mapping
quality class, read length, mean base quality along the read, mismatch / insertion / deletion rates
along the read, and the cumulative identity distribution. Tested against MultiQC 1.36.dev with
`--strict`, ruff, mypy, MultiQC's code checks and its module test suite.

1. Open an issue at https://github.com/MultiQC/MultiQC/issues ("New module: SAMStat") with a sentence on
   what SAMStat adds (mapping-quality split, short and long reads) and a link to the release.
2. Test data: fork https://github.com/MultiQC/test-data and add the example reports:

   ```sh
   packaging/multiqc/make_test_data.sh <path to samstat> <test-data fork>/data/modules/samstat
   ```

   Open a pull request there.
3. Fork https://github.com/MultiQC/MultiQC, apply the module with `git am 0001-New-module-SAMStat.patch`,
   push, and open a pull request that mentions the issue and the test-data pull request. A screenshot of
   the report helps the reviewers.

Until the module is in a MultiQC release, `samstat --multiqc` (custom content) keeps working. Users
should use one or the other: with both, SAMStat appears twice in the report.
