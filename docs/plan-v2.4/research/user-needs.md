# User needs: what researchers want from read / alignment QC

Date: 2026-10-14. Method: web searches and issue-tracker reads. Search snippets are secondary; items I could not confirm are marked UNVERIFIED. The PubMed search for 2023-2026 reviews returned only one hit (PMID 39211184, not read), so no review-literature claims are made.

## 1. SAMStat's own issue tracker (4 open issues, 0 closed)
Source: https://github.com/TimoLassmann/samstat/issues
- #12 (Jan 2023, osilander): long-read user. Wants error plots as FRACTION of reads (normalised per bin) rather than counts, and a constant, user-set bin size. Default steps (1,10,50,100) are too small at ~10 kb mean read length. Still open, no maintainer reply. https://github.com/TimoLassmann/samstat/issues/12 . The reporter also says the tool is "great and easy to use".
- #13 (Mar 2023): build failure on CentOS 7 (CMake 3.9, submodule update fails); a second user hit a link error (`tlfilename`), reportedly fixed Aug 2023. https://github.com/TimoLassmann/samstat/issues/13
- #6 (Nov 2021): request for Conda availability. (2.3.0 now ships a draft bioconda recipe.)
- #1: README image (trivial).
Reading: the only substantive user request is long-read normalisation and binning; the only other pain is installation.

## 2. FastQC (the benchmark)
- 11 open issues; themes are features and report correctness, not long reads: downsampling (#160), separate R1/R2 for paired BAM (#66), CRAM input (#54), naive averaging of Phred scores (#110), empty plots in HTML (#44), wrong adapter list (#52). https://github.com/s-andrews/FastQC/issues
- Long reads: a Galaxy post found FastQC's hard-coded JVM memory cap (-Xmx250m) made very long read jobs slow; raising it sped up all tests. https://galaxyproject.org/news/2019-04-14-fastqc-long-reads
- SeqAnswers: crashes / OOM on nanopore data; `--nogroup` produces huge plots on long reads. https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/3989-fastqc-a-quality-control-application-for-fastq-data?p=280579
- Length-distribution bins on nanopore data show odd values ending in 999 in MultiQC (root cause tracked upstream, FastQC#203; MultiQC issue open June 2026). https://github.com/MultiQC/MultiQC/issues/3591
- Third-party claim (secondary, UNVERIFIED as FastQC's own statement): per-base modules assume uniform length and are unsuitable for ONT/PacBio. https://lib.rs/crates/kira-qc
Lesson: the fixed plot set is why people use it; its failure mode is long reads and big inputs.

## 3. NanoPlot (the long-read default)
https://github.com/wdecoster/NanoPlot/issues (390 issues; first page of 25 inspected)
- Dominated by crashes and dependency breakage: Plotly/Kaleido version incompatibilities (#442, #417), blank PNG/HTML (#437, #411), BrokenProcessPool (#430), KeyError on FASTQ record (#436), empty/corrupt FASTQ handling (#416), conda install problems (#422). About 17 of 25 are errors.
- Features requested: single report with all plots (#433, closed: already existed as NanoPlot-report.html, i.e. discoverability), quality along reads (#410), low plot resolution (#424).
- #433 also asks for "what to expect" guidance for new users; maintainer calls it valuable but hard because expectations differ by assay. https://github.com/wdecoster/NanoPlot/issues/433
- Slowness: ~1 GB BAM ran >2 h; KDE plots slow; `--plots`, `--downsample`, `--huge` mitigate. https://biostars.org/p/300901
Lesson: Python plotting-stack fragility is the top annoyance; a self-contained C binary avoids it. Users want one report file and interpretation help.

## 4. Qualimap bamqc
- OOM on ~90x WGS (needs `--java-mem-size`); a 150 GB BAM did not finish in 3 days on 48 cores / 500 GB, user split per chromosome. https://biostars.org/p/69666 ; https://code.rc.uab.edu/center-for-computational-genomics-and-data-science/public/quac/-/issues/48
- Google group: 2022-23 NullPointerException on v2.2.2; wrong read counts caused by secondary alignments. https://groups.google.com/g/qualimap
- Counting secondary/supplementary correctly is a recurring correctness complaint (SAMStat 2.3 already has tests for this).

## 5. cramino (new long-read BAM/CRAM QC, shows what long-read users get fast)
https://github.com/wdecoster/cramino
- Reports: read/alignment counts, yield, N50/N75/median/mean length (base-weighted N50), identity median/mean/mode (gap-compressed, from `de` tag or NM+CIGAR), optional histograms; text/JSON/TSV.
- Claim: 140 Gbase BAM in 12 min, <1 GB RAM; reads stdin and remote files.
- Secondary alignments are never counted.
Lesson: the long-read summary users expect is N50 + yield + identity + length histogram, computed streaming and fast.

## 6. MultiQC (where reports end up)
https://github.com/MultiQC/MultiQC/issues (page 1 of 9 inspected)
- log_filesize_limit 50 MB too low for high-coverage data (#3534); plots from many samples (#3606, #3662 FastQC plot issues); interest in ECharts replacing Plotly (#3535); new cramino module requested (#3603).
- Implication: per-sample outputs must be small (SAMStat's 100 MB JSON on 1 Mb reads would hit such limits) and MultiQC-friendly.

## 7. Alignment-level questions on Biostars
- Per-read identity from a BAM is a repeated question; answers are ad-hoc scripts using NM tag, bedtools or htsjdk. https://biostars.org/p/305174
- Picard CollectAlignmentSummaryMetrics gives averages only; asker wants a per-read distribution for long-read CRAM. (Search snippet; thread id not confirmed.) 
- Rule of thumb quoted: raw ONT Q~8-12, Illumina Q~30-35; NanoPlot author suggests mean Q>12 cut-off. https://www.biostars.org/p/381829 (snippet only)
- Tool landscape: chopper filters but writes no report; seqkit stats gives N50 summaries. https://seqtec.readthedocs.io/en/latest/tools/quality-control/chopper.html

## 8. Synthesis

Top needs (ranked by evidence weight)
1. Works on long reads without breaking: bounded memory and output size regardless of read length; length-aware binning (SAMStat #12, FastQC long-read threads).
2. Read-length summary: N50, yield, log-scaled length histogram (cramino, NanoPlot). SAMStat lacks all three.
3. Per-read identity / accuracy distribution from alignments (Biostars, cramino). SAMStat has mismatch/indel profiles but no per-read identity.
4. Speed on big files, bounded RAM, stdin streaming (Qualimap, NanoPlot, cramino claims).
5. One self-contained report file; fixed plot set; short "what to expect / what is bad" guidance (NanoPlot #433).
6. Normalised axes (fraction of reads) so samples with different depth compare (SAMStat #12).
7. Clean MultiQC integration with small per-sample outputs (MultiQC #3534, #3603).
8. Easy install: conda/bioconda, static binary (SAMStat #6, #13; NanoPlot #422).
9. Correct counting: secondary/supplementary, R1/R2 separation for paired BAM (Qualimap, FastQC #66).
10. Downsampling option (FastQC #160, NanoPlot --downsample); SAMStat already has -p.

Top annoyances
- Crashes from dependency / runtime stacks (NanoPlot Python plotting; FastQC/Qualimap JVM heap).
- Timeouts and OOM on large BAMs (Qualimap, NanoPlot KDE).
- Plots that are meaningless or enormous for long reads (FastQC per-position, SAMStat per-position past a few kb).
- Hidden features / not knowing which file to open (NanoPlot #433).
- Silent miscounting (secondary alignments, naive Phred averaging FastQC #110).
- Build friction on old clusters (SAMStat #13).

Not found / UNVERIFIED
- No complaint evidence found for or against Plotly-CDN (offline report) dependency; this is a plausible need for HPC nodes without internet but UNVERIFIED by user reports here.
- No user demand found for nf-core module or native MultiQC module beyond MultiQC module requests for other tools; UNVERIFIED for SAMStat.
- Biostars/SEQanswers threads were read via search snippets, not full text; counts and quotes are therefore indicative.
- Recent (2023-2026) review papers on sequencing QC were not located.

## 9. Implications for the maintainer's "simple, standard plots" principle
- The evidence supports a small fixed set; the demand is for those plots to stay valid for long reads (scaled axes, fraction-of-reads, log length), not for more plot types.
- Cheapest high-value additions: N50/yield/length histogram, per-read identity histogram, fraction-normalised error profiles with length-adaptive bins, bounded output size.
- Keep SAMStat's differentiator (error profiles by MAPQ) visible; no competitor above provides it in the default report.
