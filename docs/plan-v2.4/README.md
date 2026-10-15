# SAMStat 2.4 master plan

Date: 2026-10-14. Branch: `planning/v2.4`. Base: SAMStat 2.3.0 (merged to main, commit c139e9d, not yet tagged).

## Files

| File | Purpose |
|---|---|
| `README.md` | This overview |
| `standard-plots.md` | Binding spec: report layout, the six plots, captions, data model, JSON v2, MultiQC |
| `work-packages.md` | 23 packages in 6 waves (0-5) for the multi-agent implementation run, with file ownership, tests and gates |
| `decisions.md` | Every open decision with options and the ruling the plan assumes; sign off before gate G0 |
| `simplicity-review.md` | Cross-topic review; its cuts are binding |
| `review-log.md` | Verdict on every issue from the plan review and what changed |

Precedence: for exact text, data definitions and rendering rules `standard-plots.md` wins; for scope, budget and cuts `simplicity-review.md` wins; `decisions.md` rulings win over both until the maintainer changes them.
| `topics/*.md` | Nine topic designs (background; superseded where they conflict with the files above) |
| `research/*.md` | Research notes behind the topics |

## Guiding principle

Make SAMStat the most useful and fastest QC tool for long and short reads, and keep it simple. FastQC succeeded with a standard set of plots that are learned once and read at a glance. SAMStat 2.4 does the same: **one page, six plots, always in the same order, the same encodings for every input, no options that change the report.** Every element answers one of four questions: did the run work, short or long reads, how well did reads map, and where along the read do errors occur by mapping quality. Anything else is cut or left to other tools.

## Goals for 2.4

1. **Long reads work.** Memory and output no longer grow with read length (today one 1 Mb read costs 3.8 GB RSS and a 96 MB JSON). N50, log-scale lengths and per-read identity are added.
2. **The report is readable in 30 seconds** and works offline: inline SVG, no Plotly CDN, no JavaScript.
3. **Fast on one big file:** a reader thread plus workers per file and htslib decompression threads, all within `-t` total threads; byte-identical output for any `-t`. Speed targets are ratios measured on x86-64 Linux (`work-packages.md` WP2.3).
4. **Zero required options:** content-based format detection, `-` for stdin, CRAM via `REF_PATH`/`REF_CACHE` (dropped first if the schedule slips), fewer documented options (12 to 8).
5. **One stable contract** for pipelines: JSON `format_version` 2, frozen for 2.x; MultiQC custom content gains N50, mean Q, median identity and the error profile (on the % of read length axis for every sample).
6. **Installable:** bioconda (2.3.0 submitted now; the bot bumps to 2.4.0).

## The standard report at a glance

Header line (incl. layout and median length) · summary strip (reads, bases, N50, mean Q, % mapped, % MAPQ >= 30, % bases clipped, median identity; the last four on aligned input only) · factual notes · MAPQ table (doubles as the legend) · plots 1-4 on every input, plots 5-6 on aligned input:

| # | Title | What it shows | Input | FastQC equivalent |
|---|---|---|---|---|
| 1 | How long are the reads? | Length histogram (linear x for short reads, log x for long), N50 line | all | Sequence Length Distribution |
| 2 | How good is each read? | Per-read mean Q histogram | all | Per Sequence Quality Scores (different averaging) |
| 3 | Where along the read does quality drop? | Mean Q per position | all | Per Base Sequence Quality (different averaging) |
| 4 | Is base composition even along the read? | % A/C/G/T/N per position | all | Per Base Sequence Content + N Content |
| 5 | Where along the read do errors occur? | Mismatch, insertion, deletion rate per aligned base (3 panels) | aligned | none: SAMStat's differentiator |
| 6 | How accurate are the aligned reads? | Per-read gap-compressed identity, median line | aligned | none |

On aligned input, plots 1-3, 5 and 6 have one line per MAPQ class (>= 30, 10-29, < 10, unmapped; unmapped is not drawn in 5-6), each with a fixed colour and dash pattern; plot 4 pools all reads. Input counts as aligned when the header has `@SQ` lines or any read is mapped, so an unaligned Dorado/Revio BAM gets the unaligned layout. Qualities are averaged as error probabilities, so mean Q reads lower than FastQC's average of Q scores; the captions say so. Paired input shows R1 and R2 as side-by-side panels. The position axis is in bp for short reads and in % of read length when the median read is >= 1 kb, chosen automatically. No pass/fail flags in 2.4.

## Roadmap

| Release | Content |
|---|---|
| **2.3.0** (now) | Tag as it is after the offline tarball check; submit to bioconda. No hotfixes for long-read memory or threading. |
| **2.4.0** | Everything in `work-packages.md`: htslib bump, htslib reader for all formats (stdin, content detection; protein FASTA unsupported), bounded collector with long-read statistics, intra-file threads, inline-SVG report with six plots, JSON v2, MultiQC additions, CRAM, deprecated options as hidden no-ops. About 35 developer-days by the topic estimates; 6 waves (0-5), 23 packages. |
| **After 2.4.0** | bio.tools entry; nf-core module (after bioconda 2.4.0); MultiQC module-request issue. |
| **2.5** | Checks (OK/CHECK) and numeric caption ranges, only after calibration on real Dorado, Revio and NovaSeq data; `--reference` and libdeflate if measured or requested; CRAM if it slipped from 2.4. |
| **3.0** | Remove `-l`, `--seed`, `--norm`, `--plotend`; reconsider the default output location. |

## Non-goals

- No GC, duplication, overrepresented-sequence, adapter, k-mer or per-tile modules; no 2D plots, coverage, insert size or homopolymer analysis; no modified-base analysis in 2.4 (revisit only on user request, `decisions.md` A11).
- No interactive features (zoom, dropdowns, toggles, dark mode), no multi-sample view (MultiQC does that), no configurable MAPQ classes, bins or thresholds.
- No pass/warn/fail verdicts in 2.4.
- No new CLI options; no `--long`/`--mode`, `--reference`, `-o` rename, `--quiet` or core auto-detection.
- No native MultiQC module, Galaxy wrapper, Docker image or Windows build.
- No perf CI gates; one manual `tests/perf/bench.sh`.

## How to run the implementation

1. Maintainer signs off `decisions.md` and does gate G0 (tag 2.3.0, bioconda PR, push the tld strbuf fix).
2. A workflow runs waves 1-5 from `work-packages.md`: per package one implementer, plus one Opus reviewer for M/L packages and Opus-implemented ones (42 agent runs in total); per gate one verifier that merges in the given order and keeps `ctest` green in Release and sanitizer builds after every merge.
3. Goldens change only in the two labelled regeneration commits (gates G3 and G4).
4. Gate G5 is manual: real-data checks, browser and print check, benchmark numbers, tag 2.4.0.
