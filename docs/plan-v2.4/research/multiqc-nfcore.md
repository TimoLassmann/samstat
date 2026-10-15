# Ecosystem integration: MultiQC, nf-core, bioconda, Galaxy

Researched 2026-10-14. Sources are GitHub/docs pages fetched that day; anything not fetched is marked UNVERIFIED.

## Bottom line (priority order)

1. **Bioconda first.** It is the gate for everything else (nf-core, Galaxy, biocontainers all consume it). Checked: there is no `samstat` recipe in bioconda-recipes (`gh api repos/bioconda/bioconda-recipes/contents/recipes` has no samstat; api.anaconda.org returns "samstat could not be found"). The draft in `packaging/bioconda/meta.yaml` is the right start. Effort: about 0.5 day plus review wait.
2. **Keep the 2.3 `--multiqc` custom-content output as the MultiQC story for 2.4.** It works today with no upstream approval. Native module only after bioconda downloads and users exist (see below). Effort for native module: 2-3 days of Python, plus an unbounded review wait.
3. **nf-core module** after bioconda: about 1 day. Low risk, high reach.
4. **Galaxy: skip** unless a Galaxy admin asks. Low value for the cost.

## 1. Bioconda / biocontainers

- Recipe lives in `bioconda/bioconda-recipes/recipes/samstat/meta.yaml` (+ `build.sh`). Submitted as a PR; bioconda-bot builds linux-64, osx-64 and osx-arm64 and a biocontainer (Docker + Singularity image on quay.io/biocontainers) appears automatically after merge. Bioconda docs: https://bioconda.github.io/contributor/workflow.html (not fetched; process from general knowledge, UNVERIFIED in detail).
- Existing draft needs: `build.sh` (cmake configure/build/install; our htslib is built by ExternalProject, so conda must not fetch from the network at build time -> the bundled htslib tarball in `thirdparty/` must be inside the release tarball). Better long term: depend on conda `htslib` and skip ExternalProject, but that changes the build; do NOT do it in 2.4 unless the bundled build fails in bioconda CI.
- Source tarball must contain the `src/tld` submodule contents (GitHub auto-tarballs omit submodules). The draft already points at a custom release asset `samstat-<version>.tar.gz`; keep that.
- Add `test: commands: samstat -v` (present) and one real run on a tiny file shipped in the tarball (`samstat --json --multiqc`), so the container is proven functional.
- Add `extra: identifiers: biotools:samstat` (UNVERIFIED that a bio.tools entry exists; check https://bio.tools/samstat; register if missing, 10 minutes).
- Steps: tag v2.4.0 -> release workflow attaches tarball + sha256 -> fill sha256 -> PR to bioconda-recipes -> fix lint -> merge -> wait for `quay.io/biocontainers/samstat:2.4.0--<hash>`.
- Gotcha: the tool name `samstat` collides with nothing in bioconda (checked).

## 2. MultiQC

### Current state (2.3.0)
`--multiqc` writes one custom-content JSON per section (general stats, MAPQ table, length and base-quality linegraphs). Works with MultiQC 1.35 (latest release at research time: v1.35, https://github.com/MultiQC/MultiQC/releases). Limits: not shown as a named module in the report header, no software version row, no per-position error-profile plots, user must set `fn_clean_sample_names: false`.

### Native module: requirements (https://docs.seqera.io/multiqc/development/modules)
- Files: `multiqc/modules/samstat/{__init__.py,samstat.py,tests/test_samstat.py}`.
- Registration: entry point in MultiQC `pyproject.toml` (`[project.entry-points."multiqc.modules.v1"]`), name in `module_order` of `config_defaults.yaml`, search pattern in `multiqc/search_patterns.yaml` (key prefix = module name, e.g. `samstat: {fn: "*.samstat.json", contents: '"samstat_version"', num_lines: 20}`).
- Class extends `BaseMultiqcModule`; required `name, anchor, href, info, doi, license, license_url` (CI linter requires license). Our DOI: 10.1093/bioinformatics/btad019 (from the paper citation in the task brief; verify before use). License GPL-3.0-or-later per `packaging/bioconda/meta.yaml`.
- Must raise `ModuleNoSamplesFound`, call `add_software_version`, `general_stats_addcols`, `add_section`, `write_data_file`.
- Tests: one blanket test runs every module on the `MultiQC/test-data` repo (a PR to that repo with example files is required: `data/modules/samstat/`). CI runs `multiqc --strict test-data`, mypy, ruff.
- Policy that matters for us (https://github.com/MultiQC/MultiQC/blob/main/.github/CONTRIBUTING.md): "MultiQC modules parse log files, they don't calculate new metrics"; modules must be fast and scale to 5000 samples. So the module must parse our existing `--json`, and should summarise (a few plots) not replicate every plot. Per-position, MAPQ-stratified error profiles are exactly our value, but with 5000 samples a MAPQ-stratified per-position plot is too heavy; show only the all-reads (or MAPQ>=30) curve per sample, with the per-sample HTML linked.
- Prioritisation: new module requests are scored automatically (https://github.com/MultiQC/MultiQC/blob/main/.claude/docs/module-triage-system.md): tool popularity 25, package downloads (PyPI/Conda) 15, community engagement 35, example files 20, feasibility 15. A tool with no conda package scores 0 on downloads, and upstream recently has modules for far newer tools, so being in bioconda first is a prerequisite in practice.
- Note: a module author is welcome to submit a PR directly without waiting for triage (the contributing guide describes PRs generally; UNVERIFIED whether a PR without prior issue is accepted for a tool with few stars). Open an issue with example JSONs first; it costs 15 minutes.

### Effort and recommendation
- Python module: ~150-250 lines reading our JSON (schema in `docs/samstat.schema.json`). Plots: general stats (reads, %MAPQ>=30, %unmapped, %mapped-identity if added), MAPQ stacked bargraph, read-length linegraph (log-binned), mean base quality per position, mismatch/ins/del rate per position. 2-3 days including test-data PR and strict-mode fixes; review latency is the unknown (weeks to months, UNVERIFIED).
- Do it in the order: (a) release 2.4 in bioconda, (b) open MultiQC module-request issue with example JSON, (c) write the module only if the issue gets traction, or the maintainer is willing to carry a Python module. Until then, custom content stays. This is also the lowest-maintenance path: the custom-content writer (`mqc_report.c`) already exists and is tested.
- Simplification option: since a native module will read the same JSON, once it lands consider deleting `mqc_report.c` (-1 output format to maintain). Do not delete before the module ships in a MultiQC release, and keep it for a deprecation cycle.
- Design consequence for 2.4: make the JSON the single machine-readable contract (small for long reads: pre-binned, bounded size). Every integration (MultiQC module, nf-core, anything else) then reads one stable file. Add `samstat_version` and a `schema_version` key so a MultiQC search pattern and version reporting are trivial (check that JSON already has them; UNVERIFIED, not inspected).

## 3. nf-core module

Spec: https://nf-co.re/docs/guidelines/components/modules ; template inspected: `modules/nf-core/fastqc` and `modules/nf-core/cramino` (long-read BAM stats tool, a close analogue) in https://github.com/nf-core/modules. Checked: no `samstat` module exists. Other QC modules exist (fastqc, fastp, fastplong, nanoplot, chopper, cramino, seqkit, samtools), so SAMStat competes in a crowded list; justify by the MAPQ-stratified error profile.

Module layout `modules/nf-core/samstat/`:
- `main.nf`: process `SAMSTAT`, `tag "${meta.id}"`, `label 'process_low'` (or `process_medium` with `-t`), `conda "${moduleDir}/environment.yml"`, `container` directive (Seqera Wave community container for docker and singularity; generated by the nf-core tooling, no hand-written URL).
- Input `tuple val(meta), path(reads)` (BAM/SAM/FASTQ/FASTA). Outputs: `html`, `json`, `mqc` (the `_mqc.json` files, glob `*_mqc.json`), and versions.
- Versions: current convention is the topic channel, not `versions.yml`: `tuple val("${task.process}"), val('samstat'), eval('samstat -v | sed ...'), emit: versions_samstat, topic: versions`. Version string must start with a digit (no `v`). So `samstat -v` output must be parseable by sed; make the output a clean `samstat 2.4.0` or just the number. Concrete 2.4 requirement: stable `samstat -v` format.
- `when: task.ext.when == null || task.ext.when`; options via `task.ext.args`; `prefix = task.ext.prefix ?: "${meta.id}"`.
- `stub:` block required: `touch ${prefix}.html`, `touch ${prefix}.json`, `touch ${prefix}.samstat_general_mqc.json` etc.
- `meta.yml`: name, description, keywords, tools (with homepage, documentation, licence ["GPL-3.0-or-later"], identifier `biotools:samstat`), `input`, `output` (each channel with ontologies, pattern), `topics: versions`, `authors`. See fastqc `meta.yml` for the exact shape (inspected).
- `environment.yml`: `channels: [conda-forge, bioconda]`, `dependencies: [bioconda::samstat=2.4.0]`.
- `tests/main.nf.test` + `main.nf.test.snap` (+ `nextflow.config` for ext.args): test with `params.modules_testdata_base_path` files from nf-core/test-datasets (existing small human BAM/FASTQ there are reusable, so no new test data needed, UNVERIFIED which); a `-stub` test; snapshot only deterministic outputs. Our outputs have no timestamps/absolute paths (stated in `docs/multiqc.md`), so snapshot md5s are stable: a nice side effect of the 2.3 work. HTML includes version only; fine to snapshot a name, not content.
- Process: `nf-core modules create samstat` scaffolds everything; `nf-core modules lint samstat` and `nf-core modules test samstat` (nf-test, docker/conda profiles) locally; PR to nf-core/modules; 2 reviews typically from the community (Slack #modules for questions). Nothing here is blocked on MultiQC native support.
- Effort: ~1 day, once bioconda has 2.4.0 (the module test pulls the conda package/container). Pipelines then wire `SAMSTAT.out.mqc` (or `json`) into the MultiQC process: with custom content that works today, no change needed in MultiQC.

## 4. Galaxy

- Wrapper in https://github.com/galaxyproject/tools-iuc (`tools/samstat/` with `samstat.xml`, `.shed.yml`, `test-data/`, planemo tests). Tools use a conda requirement `<requirement type="package" version="2.4.0">samstat</requirement>` and are auto-tested by planemo in CI; on merge uploaded to the Tool Shed (https://github.com/galaxyproject/tools-iuc). No samstat tool in the current IUC listing I could see (fastqc, nanoplot, cramino exist; the listing page did not enumerate tools, so absence is partly UNVERIFIED; `gh api .../tools` grep for samstat returned nothing).
- Effort: ~1 day XML + tests; ongoing: each release needs a version bump PR (IUC bots help). Recommendation: do not write it; accept a community PR if one appears. A Galaxy-side maintenance cost conflicts with the single-maintainer constraint.

## 5. Things 2.4 must do so integrations stay cheap

1. `samstat -v` prints a bare parseable version (nf-core versions topic, bioconda test, MultiQC).
2. Release tarball includes submodule and bundled htslib (bioconda, offline build).
3. JSON is the stable contract: versioned, bounded size (long reads pre-binned), deterministic. The HTML and MultiQC outputs are views on it.
4. Exit codes 0/1/2 already documented; nf-core handles non-zero as failure, fine.
5. `-t` maps to `task.cpus`: nf-core module passes `-t ${task.cpus}`; note that today a single BAM uses only one core, so `label 'process_single'` is honest until the single-file parallelism work (other research doc) lands.
6. No network at runtime: Plotly from CDN breaks offline HPC nodes and containers without internet; embed Plotly (or lighter) in the HTML. nf-core pipelines run on offline clusters, so this is a real integration issue.
7. Output file names deterministic from the input name and `-d`; the nf-core `prefix` handling needs predictable names (or an explicit output-prefix option). Check current naming (`<file>.html` etc.; UNVERIFIED how names behave for `.bam` vs staged symlinks).

## 6. Suggested sequence and cost

| step | effort | depends on | value |
|---|---|---|---|
| tag 2.4.0, tarball, sha256 | 0.5 d | code done | gate |
| bioconda PR | 0.5 d + wait | tarball | high (container, conda) |
| bio.tools entry | 0.1 d | none | low |
| nf-core module | 1 d | bioconda | high |
| MultiQC module request issue + example JSON | 0.2 d | JSON stable | medium |
| MultiQC native module PR | 2-3 d + review | issue traction, bioconda downloads | medium (nicer plots; custom content already works) |
| Galaxy wrapper | 1 d | bioconda | low, skip |
