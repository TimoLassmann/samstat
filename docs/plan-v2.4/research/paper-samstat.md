# SAMStat 2 paper (Lassmann 2023): claims, plots, what 2.4 should strengthen

Source: Lassmann T. "SAMStat 2: quality control for next generation sequencing data." Bioinformatics 39(1), btad019 (2023). PMID 36637208.
URL: https://academic.oup.com/bioinformatics/article/39/1/btad019/6986976 (Applications Note; read via page fetch, not the PDF, so wording is from a summarised page: quotes below should be re-checked against the PDF before reuse in a 2.4 paper/changelog - UNVERIFIED verbatim).
Local sources: Readme.org, CHANGELOG.md, docs/golden-diff-v2.3.md.

## 1. Unique value claimed by the paper

- Core claim (quote): SAMStat "is unique by focusing on revealing differences between well and poorly mapped reads." Abstract: it "displays sequence composition, base quality composition and mapping error profiles split by mapping quality."
- Rationale: flagstat-type numbers show mapping is poor but not why (adapters, N bases, low quality at positions, aligner settings). MAPQ is used as a phred-scaled probability that the alignment is wrong, so contrasting high-MAPQ vs low-MAPQ vs unmapped reads exposes the cause.
- Positioning: "complements existing tools" (FastQC, MultiQC, RNA-SeQC 2 are cited; no benchmark or head-to-head comparison).
- Not claimed as unique: speed, format support. Those are 2.0 improvements.

## 2. What changed in v2 (as described)

- C rewrite; claimed "at least twice as fast" (no benchmark shown).
- htslib API instead of samtools+pipes; FSM FASTA/FASTQ parser.
- Paired-end (per-fragment split) and long-read support (long reads/names).
- Multi-file parallelism via threads.
- Interactive Plotly plots (zoom, motivated by long reads); read-end views and full-length views.
- Single HTML report.

## 3. Plots highlighted

Report contents: basic statistics; base composition (split by fragment); length distribution; mismatch, insertion, deletion profiles along read; buttons to toggle MAPQ bins. Soft clips are not shown.

Only one figure (3 panels), all demonstrating MAPQ stratification on ENCODE BAMs (ENCFF039WAK, ENCFF050BNR, ENCFF596PVI):
- (a) base composition, polyA+ RNA-seq: N spike at read position 40.
- (b) length distribution, long-read RNA-seq: unmapped reads bimodal (~165 nt and ~1360 nt), mapped mostly >250 nt; MAPQ bins differ in length profile.
- (c) mean base quality, 10x single cell: unmapped reads clearly lower quality than aligned.

Observation: the paper's three showcase diagnoses are base composition, length and base quality, each split by MAPQ. The mismatch/insertion/deletion profiles are listed but not shown in the figure, although the Readme calls them out as a headline output.

## 4. Weaknesses in the paper (what a reviewer or user would hit)

- Speed claim ("2x faster") has no benchmark. Single BAM still uses one core in 2.3.0.
- No accuracy/validation; only 3 ENCODE files, and only qualitative reading.
- Internal inconsistency: text says Nanopore MinION, Fig 1b caption says PacBio (UNVERIFIED which is right).
- Long-read support is claimed but 2.3.0 has known long-read weaknesses (1 Mb read -> ~100 MB JSON, per-position plots meaningless past a few kb, no N50/log-length histogram/identity).
- Soft clips not displayed (adapter diagnosis is the headline use case, and adapters often sit in soft clips; 2.3.0 now excludes them from per-position stats, per CHANGELOG).
- MAPQ dependence: the method's validity rests on aligner MAPQ semantics (minimap2, bwa, STAR differ; some aligners give few distinct values). Four fixed bins (>=30, 10-29, <10, unmapped) are not justified in the paper.
- Offline use: Plotly from CDN; report dies without internet (HPC nodes).
- 2.2 numerical errors since fixed in 2.3.0 (CHANGELOG): FASTQ quality offset 0 instead of 33 (every quality +33), soft-clip handling at 5' end, supplementary alignments counted as reads, buffer overflows. Any paper figure produced with 2.2 FASTQ quality would have been wrong; the claim "unmapped reads have lower quality" came from BAM input so likely unaffected, but the position-40 N spike and length plots should be re-run on 2.3 (UNVERIFIED whether numbers change).

## 5. Claims 2.4 should strengthen (ranked, simplicity-first)

1. MAPQ-stratified view stays the signature. Make it the default, not an option: same stratified plots every run, same order, with a one-line "how to read this" per plot. Keep the 4 bins; do not add more.
2. "Why is mapping poor" in one glance: a top-of-report summary showing the unmapped / low-MAPQ vs high-MAPQ contrast for the three paper-proven signals (composition, length, quality) plus error profiles. Reuse existing data; no new metrics beyond what makes long reads interpretable.
3. Long reads actually work: log-scaled length histogram + N50, error profiles on binned or capped positions (e.g. first/last few kb), bounded JSON/HTML size. This turns the "supports long reads" claim into a true one. Priority over any new feature.
4. Substantiate speed: add a reproducible benchmark (time, RSS) vs FastQC/samtools stats on short and long files, and parallelism within a file (the paper only has per-file threading). Do not claim "2x" again without data.
5. Reproducibility: ship the paper's three ENCODE examples as a documented regression/demo set (golden reports) so every claim in the paper is reproducible by one command (the paper says reproduction code is on GitHub; confirm it still matches 2.3 - UNVERIFIED).
6. Offline-capable report (vendor Plotly or lightweight static plots), so the report is a single self-contained file, as FastQC's is.
7. Fit into pipelines without adding scope: JSON + MultiQC already exist in 2.3.0 and extend the "complements FastQC/MultiQC" claim; a native MultiQC module would make the complement concrete (check scope against maintainer budget).
8. Resolve soft-clip story: either show clipped-base composition/length (adapter evidence) or state clearly why not. Candidate small addition: soft-clip length distribution per MAPQ bin.

## 6. Things to avoid (per guiding principle)

- No new plot families beyond the standard set; no per-aligner MAPQ calibration options; no configurable bin schemes.
- Do not cite "2x faster" or "long-read support" in 2.4 material until backed by benchmark / fixed behaviour.

## 7. Items not verified

- Verbatim paper text (fetched via summarising tool, PDF not read).
- Nanopore vs PacBio for ENCFF050BNR.
- Whether paper reproduction code still runs against 2.3.
- License line (not stated on the fetched page).
