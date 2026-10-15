# Short-read QC in 2026: what SAMStat 2.4 should (and should not) do

Scope: fastp, Falco, MultiQC conventions, FastQC module value, binned Q-scores on modern platforms. Evidence marked UNVERIFIED where not confirmed from a source.

## 1. Findings from the landscape

### fastp (https://github.com/OpenGene/fastp)
- An all-in-one preprocessor, not primarily a QC viewer: JSON + HTML report, before/after filtering comparison (quality curves, base content, k-mer, Q20/Q30, GC, duplication, adapter content).
- Duplication: on by default, ~1 GB RAM, +10-20% runtime. Duplicate = all bases identical (errors/Ns hide duplicates). Hash method may misclassify ~0.01%.
- Overrepresented sequence analysis is OFF by default (`-p`), sampled 1/20 reads. The maintainers themselves treat it as optional.
- Adapters: auto-detected for single-end; paired-end uses overlap analysis (auto-detect optional). PolyG trimming automatic for NovaSeq/NextSeq (detected via machine ID); polyX opt-in.
- Long reads: README directs to a separate tool, fastplong (the README is internally inconsistent about long-read support). Version not confirmed from the page (mentions v1.4.0 in an example and "fastp 1.0" in citation) - UNVERIFIED.
- Lesson: users value one fast command producing one JSON + one HTML with a small fixed set of numbers (Q20/Q30, GC, dup rate, adapter %, reads before/after). Speed is a headline feature.

### Falco (https://github.com/smithlabcode/falco)
- v2.0 is a FastQC drop-in re-implementation: same three outputs (`fastqc_data.txt`, `fastqc_report.html`, `summary.txt` with pass/warn/fail per module), so MultiQC reads it as FastQC.
- Deviations that show where FastQC modules are weakest: tile analysis uses ~1/10 of reads and grade can flip between runs with multithreading; duplication hashes 1,000,000 reads spread across the file instead of the first 100,000 (old method via `--orig-dups`); tile analysis disabled for stdin unless tile field given.
- Lesson: the FastQC format is the de-facto interchange standard; reproducing its text format earned Falco adoption without new plots. Sampling choices (first N reads vs spread) matter for correctness.

### MultiQC conventions
- MultiQC General Statistics is deliberately narrow ("keep it as narrow as possible" - maintainer statement, https://community.seqera.io/t/picard-qualityyieldmetrics-doesnt-display-all-data-in-qc-report/1805). Exact fastp/FastQC general-stats column lists could NOT be verified from docs fetched (https://docs.seqera.io/multiqc/modules/fastp, .../fastqc). From memory (UNVERIFIED): FastQC contributes %Dups, %GC, avg length, %Fail, M Seqs; fastp contributes %Dups, %GC, %PF, M Seqs and adapter/filter info.
- FastQC module in MultiQC shows pass/warn/fail status heatmap across samples, top overrepresented sequences table, theoretical GC curve (hg38/mm10). Cross-sample comparison is the main use of QC at scale; single-file HTML is secondary.
- SAMStat 2.3 already emits `_mqc` custom content. A thin, stable set of general-stats columns (reads, %GC, mean length, %Q30, %dup if cheap, %mapped, %MAPQ>=30) is the highest-leverage integration; per-sample plots are lower value.

### Binned and new-platform qualities
- NovaSeq X: 3 Q bins. Control Software v1.0/1.1: Q12 (3-17), Q20 (18-29), Q37 (30+). v1.2+: Q12, Q24, Q40 (https://knowledge.illumina.com/instrumentation/novaseq-x-x-plus/instrumentation-novaseq-x-x-plus-reference_material-list/000008320). One Illumina article labels the middle bin Q23 (https://www.illumina.com/science/genomics-research/articles/data-quality-q-scores.html); the bin label is the average of the range, not a threshold. Consequence: Q30 fraction computed from binned data is a bin-boundary artefact (a base labelled Q24 counts as <Q30 though its range spans 18-29; Q37/Q40 is >=30).
- Per-base quality boxplots on 3-level data degenerate to flat lines/steps; mean quality per position is still computable but the distribution is near-meaningless. I found no source saying FastQC is "uninformative" on binned data; this is an inference (UNVERIFIED as a community consensus). One Biostars user reports a normal-looking plot (https://biostars.org/p/9597335).
- Element AVITI: continuous (unbinned) integer Q-scores; legacy tools built for binned input can drop reads (https://www.elementbiosciences.com/knowledge-base/does-the-aviti-bin-quality-scores). Vendor claims Q40+; UltraQ claims 90% reads >=Q40 (https://www.elementbiosciences.com/blog/enhancing-sequencing-accuracy-with-cloudbreak-ultraq) - vendor claims. Implication: the parser must not assume binning, and quality range may exceed SAMStat's Q41 cap (check: UltraQ reports Q50+; cap at 41 would clip - flag for tld/collect review).
- Ultima UG 100: I could not verify binning or FASTQ/CRAM quality conventions. UNVERIFIED. It is CRAM-native and homopolymer-error dominated (UNVERIFIED), which means CRAM input and homopolymer-length-aware indel error matter more than per-base Q plots.
- MGI/DNBSEQ: UNVERIFIED, not researched; no data gathered.
- NovaSeq patterned flowcell duplicates: only forum anecdote found (https://www.seqanswers.com/forum/sequencing-technologies-companies/illumina-solexa/60394-novaseq-from-illumina); optical vs PCR dup separation needs read-name tile/x/y parsing. No authoritative guidance found - UNVERIFIED.

## 2. FastQC module value assessment in 2026

| FastQC module | Still useful? | Note for SAMStat |
|---|---|---|
| Basic statistics | Yes | Keep: reads, length, GC, %Q30 (caveat on binned) |
| Per-base quality | Declining on binned platforms | SAMStat already shows quality by position; stratify by MAPQ is the differentiator. Show mean + Q-bin fractions rather than boxplot |
| Per-tile quality | Niche, Illumina-only, flaky (Falco sampling) | Skip. Needs read-name parsing, platform-specific |
| Per-sequence quality | Weak on binned data (few modes) | Merge into one "quality distribution" bar, or drop |
| Per-base sequence content | Misleading for RNA-seq/amplicon/random-hexamer (always warns); useful for spotting adapters/primers | Keep as the one composition plot (already exists) |
| Per-sequence GC | Useful: contamination, bias | Keep |
| Per-base N | Low value; one line in basic stats | Fold into per-base content (N as fifth line) or report N% only |
| Length distribution | Essential, and the foundation for long reads | Keep, log scale; same plot serves short and long |
| Duplication levels | Useful but expensive and sampling-sensitive; first-100k method biased | Optional/cheap estimate only (sampled hash), report one number. UNVERIFIED value claim |
| Overrepresented sequences | Useful for adapters/contamination, but fastp keeps it opt-in | Fold in: report top few k-mers/sequences that are adapters; not a separate heavy module |
| Adapter content | Highly useful for short reads | Keep: simple cumulative % of reads containing known adapter prefix per position, built-in list of ~5 adapters |
| K-mer content | Little used; noisy | Drop (already absent) |

## 3. What short-read users still need from a QC tool
1. One command, zero required options, sub-minute on a typical lane; multithread within a single file (SAMStat is single-core per file: a gap vs fastp/Falco).
2. A consistent small plot set (read length, per-position quality, per-position composition, GC, adapter content, MAPQ/mapped fractions and error profiles for alignments).
3. Machine-readable summary with stable keys for MultiQC and pipelines (SAMStat 2.3 has JSON + _mqc).
4. Pass/warn/fail at module level (FastQC/Falco `summary.txt` convention) - users and pipelines gate on it. SAMStat has none; if added, keep thresholds few and documented. (Judgement call for the maintainer: statuses are the most-copied FastQC feature but also its most-criticised; if included, make them advisory.)
5. Offline HTML: Plotly via CDN breaks on HPC compute/air-gapped nodes. fastp/Falco reports are self-contained (UNVERIFIED for fastp; Falco HTML believed self-contained). A static, bundled or minimal inline-SVG report is a strong usability win.
6. Platform-aware defaults: handle binned (NovaSeq X, 3 levels) and continuous (AVITI, possibly >Q41) quality without user flags; report Q-bin counts instead of pretending to continuity.
7. stdin/streaming input (fastp supports STDIN/STDOUT) for pipelines.

## 4. Recommendations (ordered by value/cost)
1. Adapter content plot for FASTQ (short reads) - small built-in list, no config. High value, low complexity.
2. Quality display that degrades gracefully: if distinct Q values <= ~5, render as stacked bin fractions per position; else the continuous plot. Auto-detected, no option. Raise/remove Q41 cap or at least detect saturation.
3. Report %Q30 with a binned-data note, and a short "read this" caption under every plot (interpretability principle).
4. Multithread within one file (reader thread + worker partitions) - speed is the main reason people use fastp/Falco.
5. Self-contained HTML (no CDN).
6. Stable general-stats columns for MultiQC; consider a native MultiQC module later only if upstream is willing (UNVERIFIED effort).
7. Cheap duplication estimate (sampled hash of first-N bases + length) reported as a single number; do not build FastQC's histogram.

Explicitly do NOT add: per-tile quality, k-mer content, trimming/filtering (fastp's job - SAMStat is a QC tool), UMI handling, paired-end overlap analysis, optical duplicate detection.

## 5. Open questions / UNVERIFIED
- fastp current version and 1.x notes; fastp/FastQC general-stats exact columns.
- Ultima UG 100 and MGI quality encodings.
- Whether Falco/fastp HTML is fully offline.
- Authoritative optical-duplicate guidance for NovaSeq X.
- Community consensus on FastQC per-base quality utility with binned data.
