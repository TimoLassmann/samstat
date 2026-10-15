# Can one standard plot set serve short and long reads?

Status: research note for SAMStat 2.4 planning (2026-10-14). Short answer: yes, for roughly 9 of 10 plots. The exception is position-along-read, which needs one automatic rule (normalise the x axis), not a mode switch.

## 1. What the reference tools do

| Tool | Approach to length | Source |
|---|---|---|
| FastQC | Fixed 12 modules, short-read oriented. Length distribution is a plain histogram. | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/ (12 modules confirmed). Behaviour on long reads (interval grouping of positions, memory/slowness, pass/warn/fail thresholds tuned for Illumina) is from maintainer memory: UNVERIFIED. |
| NanoPlot | Long-read specific. Plots: read-length histogram, log-transformed length histogram, length vs quality bivariate, cumulative yield, reads per channel heatmap, length/quality over time, and (alignment input) identity vs length, aligned vs sequenced length, MAPQ vs length. N50 drawn as a mark on the histogram. `--loglength` option. | https://github.com/wdecoster/NanoPlot |
| cramino | Text/JSON/TSV table, not plots: N50, N75, median, mean length, yield, coverage, identity (median, mean, mode), aligned-only length stats. Identity is gap-compressed identity from `de` tag, else NM+CIGAR. Optional length and accuracy histograms (`--hist`), `--scaled` weights bins by bases. Claims 140 Gbase BAM in 12 min, <1 GB RAM. | https://github.com/wdecoster/cramino |
| samtools stats | Text sections, no plots itself: SN summary, FFQ/LFQ per-cycle quality, GCF/GCC, IS insert size, RL/FRL length distribution, MAPQ, ID indel size, IC indel per cycle, COV, GCD. plot-bamstats turns this into figures (its exact plot list not verified here). | https://www.htslib.org/doc/samtools-stats.html |

Lessons:
1. Long-read tools converge on **length** (log axis + N50), **yield-by-length** (cumulative or base-weighted), and **identity/accuracy per read**. These are the questions long-read users ask first.
2. Short-read tools converge on **per-position** quality and base composition. Long-read tools mostly drop these; no one finds position 37,512 of a nanopore read meaningful.
3. NanoPlot grows by adding bivariate plots; this is the bells-and-whistles path to avoid.
4. cramino shows that a one-line numeric summary (N50, median identity) covers much of what researchers need, cheaply.

## 2. What SAMStat 2.3 has today (src/report/stat_report.c)

Per R1/R2: length distribution; base composition; base-quality histogram; base composition and quality at read end; base composition and quality along the (normalised) read; mismatch / insertion / deletion profiles stratified by 4 MAPQ bins. Plus file and mapping statistics tables.

Gaps for long reads: no N50, no log length axis, no yield weighting, no per-read identity. Per-position plots past a few kb are meaningless and (via JSON) enormous.

## 3. Design principle: adapt data, not plot list

The plot set is fixed. Three automatic adaptations, all derived from observed data, none user-visible switches:

- **A1. Log-spaced length bins.** Always bin length on a log scale (about 40 bins, e.g. 10 bp to max). For 150 bp Illumina reads this is a single tall bar at the right spot, which is itself informative (uniform length = trimmed/untrimmed). Cheaper alternative: linear when max/min < ~4, log otherwise. Prefer always-log for one code path, but check that fixed-length data still reads well (open question O1).
- **A2. Position axis in fixed number of bins.** Per-position plots use at most ~50 bins. Short reads: 1 bp/bin up to 50 bp, then existing FastQC-style intervals. Long reads: positions are shown as **relative position (0-100% of read)** plus start/end panels of fixed bp (first/last 50 bp). Rule: if max length > ~1000 bp, x axis is normalised. SAMStat already has "along the read" normalised plots and "at end of reads" plots; this is mostly choosing which to show.
- **A3. Weight toward bases for long reads.** Report N50 as a number and a base-weighted length histogram (cramino `--scaled`). Reads-weighted and base-weighted histograms are the same plot with two traces (toggle in Plotly legend), not two plots.

## 4. Candidate unified sets

### Candidate A: "Ten plots" (recommended)

Every report, every file, same order. Each is per file; R1/R2 are traces in the same plot rather than separate plots (halves plot count for paired data).

| # | Plot | Short-read reading | Long-read reading | Adaptation |
|---|---|---|---|---|
| 0 | Summary strip (not a plot): reads, bases, mean/median length, **N50**, mean Q, mapped %, **median identity** (if aligned) | same | same | none |
| 1 | Read length histogram (reads and bases traces) | spike = fixed length | spread, N50 marked | A1, A3 |
| 2 | Mean quality per read histogram | per-seq quality | per-read Q (Q7-Q25 typical) | x range from data; cap Q41 -> Q60 (see 5) |
| 3 | Quality along read (heatmap or mean+IQR line) | FastQC per-base quality | decay/pattern across relative position | A2 |
| 4 | Base composition along read | adapter/bias | mostly flat; shows polyA/ends | A2 |
| 5 | GC content per read histogram | contamination | contamination | none |
| 6 | Mapping quality stratified bins (MAPQ bar: >=30, 10-29, <10, unmapped) | same | same | none |
| 7 | Mismatch profile along read by MAPQ | SAMStat signature | same, normalised | A2 |
| 8 | Insertion profile along read by MAPQ | same | indels dominate in ONT/older PacBio | A2 |
| 9 | Deletion profile along read by MAPQ | same | same | A2 |
| 10 | Identity per read histogram (aligned data) | near 99-100% | 85-99%, key long-read metric | none; derived from NM or MD + CIGAR |

Notes:
- Plot 10 is the one genuinely new plot. Cost: one accumulator per read. Also gives an error-rate view that short-read users find useful.
- Plots 7-9 could merge into one "error profile" plot with three panels, keeping the count at 8 sections visually (worth a design decision; see O3).
- Dropping/merging versus 2.3: "at end of reads" composition and quality plots (4 plots per read) can fold into plots 3/4 as an end zoom, or be dropped. Recommend dropping from the default report and using A2's relative-position axis; this removes 4 plots.
- Unaligned FASTQ/FASTA: plots 1-5 only (and 10 omitted). Same layout, fewer sections. Matches FastQC expectations.

Trade-offs: simple and stable; plot 3/4 with relative position lose the "first 20 bp" resolution for short-read adapter detection unless A2's rule keeps 1 bp bins for reads <= 1000 bp (it does).

### Candidate B: "Eight plots" (minimal)

Candidate A without plots 5 (GC), 4 (composition, keep only in short mode) and with 7-9 merged into one three-panel error figure. Pros: smallest. Cons: loses GC contamination check and composition-bias diagnostic that FastQC users expect; composition is where adapters and random-hexamer bias show up.

### Candidate C: "FastQC clone plus long tail"

FastQC's modules plus N50/identity. Rejected: duplicates FastQC without SAMStat's distinguishing value (error profiles by MAPQ), and FastQC's duplication/overrepresented/kmer modules need whole-file hashing, expensive for long reads and big BAMs.

### Candidate D: NanoPlot-like bivariate set

Length-vs-quality, identity-vs-length, MAPQ-vs-length scatter/hexbin. Powerful for ONT troubleshooting but needs 2D accumulators, heavier JSON/HTML, harder to read at a glance. Reject as default; at most keep one 2D plot (length vs mean Q) as a possible later addition (O2).

## 5. Details that matter when implementing

1. **Quality cap**: samstat caps base quality at Q41. Modern PacBio HiFi reports Q up to 93 (HiFi/CCS reads commonly Q30-Q60+ per base, UNVERIFIED for exact max) and ONT per-read Q is typically below 30. Cap should be raised (e.g. 93, the SAM maximum for printable Phred+33) with binning; otherwise HiFi quality plots saturate.
2. **Memory/JSON size**: per-position accumulators must be bounded (fixed bin count, A2). That alone fixes the ~100 MB JSON for 1 Mb reads; accumulate into relative-position bins at parse time, not into per-bp arrays.
3. **Identity definition**: use gap-compressed identity as cramino does (matches ecosystem; `de` tag if present, else NM+CIGAR) to be comparable with cramino/NanoPlot output; SAMStat already parses MD, so alternatively compute from MD+CIGAR and state the definition in the report footer. Document which.
4. **Detecting regime**: no switch exposed; the only data-driven decision is A2's threshold on max read length (1000 bp). Use a single named constant. Single-pass parse needs the max length before binning; solutions: bin on relative position always and, for short reads, additionally keep 1 bp bins up to first 100 positions (cheap, fixed size). This avoids a two-pass scan.
5. **Reporting**: the summary strip with N50 and median identity is the most valuable addition per line of code; it also feeds MultiQC general-stats.
6. **Offline HTML** is orthogonal but affects presentation: log-axis, hover and legend toggles come from Plotly; keep.

## 6. Recommendation

Adopt Candidate A with: summary strip + 10 plots, R1/R2 as traces, log length bins, relative-position axes for reads > ~1 kb, per-read identity added, "at end of reads" plots dropped, Q cap raised. Net plot count goes down for paired data (from about 15 per read-pair to 10) while gaining long-read usefulness. No `--long` flag.

## 7. Open questions for the maintainer

- O1. Always-log length axis versus automatic linear/log: test on fixed-length 150 bp data for readability.
- O2. Include length-vs-quality 2D plot? Default: no.
- O3. Merge mismatch/insertion/deletion into one three-panel figure, or keep three separate plots?
- O4. Is per-read identity plot only for aligned input (assumed yes)?
- O5. Keep 4 MAPQ bins for long reads (supplementary/secondary alignments are common; confirm they are excluded or labelled).

## 8. Unverified items
- FastQC long-read behaviour (grouping, runtime, thresholds).
- HiFi maximum base quality values.
- plot-bamstats exact plot list.
- NanoPlot and cramino current version numbers (pages do not state them).
