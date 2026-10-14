# FastQC / Falco: research for SAMStat 2.4

Date: 2026-10-14. Thresholds below were read from the official FastQC help pages (base URL
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/ ) unless marked UNVERIFIED.

## 1. The twelve FastQC modules and their flags

| # | Module | Plot / content | Warn | Fail |
|---|--------|----------------|------|------|
| 1 | Basic Statistics | table: file, type, encoding, total sequences, poor-quality flagged, length, %GC | never flagged | never flagged |
| 2 | Per base sequence quality | box/whisker of Phred per position | lower quartile <10 or median <25 at any base | lower quartile <5 or median <20 |
| 3 | Per sequence quality scores | histogram of per-read mean Q | most frequent mean Q <27 (0.2% error) | <20 (1% error) |
| 4 | Per base sequence content | A/C/G/T % lines per position | A-T or G-C difference >10% at any position | >20% |
| 5 | Per sequence GC content | histogram vs fitted theoretical normal | summed deviation >15% of reads | >30% |
| 6 | Per base N content | % N per position | N >5% at any position | >20% |
| 7 | Sequence length distribution | histogram of lengths | reads not all the same length | any zero-length read |
| 8 | Duplicate sequences | % remaining after dedup by duplication level | non-unique >20% of total | >50% |
| 9 | Overrepresented sequences | table of sequences with source hit | any sequence >0.1% of total | >1% |
| 10 | Adapter content | cumulative % reads with adapter by position | any sequence in >5% of reads | >10% |
| 11 | Kmer content | enriched 7-mers by position | binomial p <0.01 | p <1e-5 |
| 12 | Per tile sequence quality | heatmap tile x position, Illumina only | tile mean Q >2 below the position mean | >5 below |

Page quirks: the "Failure" text for modules 10, 11, 12 literally says "warning" (doc typo, intended as fail).
Duplicate sequences: only sequences first seen in the first 100,000 reads are tracked, counted to end of file;
reads >75 bp are truncated to 50 bp before comparison (official page).
Also in the overview page: a summary tab with a green tick / orange ! / red cross per module.

Source pages: Basic index https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/
(modules 2 to 12 pages, same folder, e.g. "2 Per Base Sequence Quality.html").

UNVERIFIED (from memory): per-base plots group positions for reads >50 bp (first 9 individually, then bins of
growing width; `--nogroup` disables); FastQC is Java, needs a JRE; `--nano` flag exists for Nanopore fast5 input;
FastQC samples only the first 200,000 reads for the Kmer module (kmer sampled 1 in 50); FastQC reads full file
for most other modules; output is a self-contained HTML (images embedded as base64) plus `fastqc_data.txt` and
`summary.txt` in a zip.

## 2. Falco

Source: https://github.com/smithlabcode/falco (verified via fetch).
- "An emulation of the popular FastQC software"; rewritten in 2.0; multithreaded C++; same three outputs
  (`fastqc_data.txt`, `fastqc_report.html`, `summary.txt`), so MultiQC and downstream scripts keep working.
- Differences in the README: tile analysis uses ~1/10 of reads; new duplication method hashes 1,000,000 reads
  (first 50 nt) sampled across the file, vs FastQC's first 100,000 unique ("numbers differ dramatically");
  `--orig-dups` restores the old way; stdin works but needs `--stdin fq:4` style hint for tile IDs;
  multithreading makes tile and duplication grades vary slightly run to run.
- No timing numbers on the README (they point to GitHub Discussions). Speed claims from the paper (Cancer Bioinformatics /
  F1000Research 2021, de Sena Brandine and Smith) are UNVERIFIED here.
- Lesson: the market accepts a faster re-implementation ONLY if it keeps the FastQC output contract identical.
  Falco's adoption driver is "drop-in", not new plots.

## 3. Why FastQC is the de-facto standard (synthesis; opinion, UNVERIFIED as a sourced claim)

1. Same report every time: fixed order, fixed names ("Per base sequence quality"), learned once, quotable in
   methods sections and core-facility reports.
2. Traffic-light summary: a non-expert gets a verdict in 2 seconds; an expert reads the plot.
3. Zero configuration: `fastqc file.fq.gz` gives everything; GUI for bench biologists.
4. Self-contained HTML, one per file, easy to email.
5. Each module has a help page (what good looks like, what bad looks like, why).
6. Machine-readable `fastqc_data.txt`/`summary.txt` -> MultiQC became the aggregation layer, which cemented it.
7. Early and free (Babraham), bundled into pipelines (nf-core, Galaxy, core facility defaults).

## 4. Common criticisms

- Flags are Illumina-genome-centric and widely mis-read: "FAIL" on per-base content for RNA-seq (random hexamer
  priming), Nextera/ATAC, amplicons, bisulfite; duplication FAIL for RNA-seq/high-depth; overrepresented and
  length warnings for trimmed data. Experienced users say "ignore the flags, read the plots".
- Thresholds are fixed, not tunable in the report, and arbitrary (limits file exists: `--limits`; UNVERIFIED detail).
- Duplication estimated from the first 100k reads; kmer module slow and rarely interpretable.
- Per sequence GC compares against a fitted normal, not the organism's genome; contaminated or high-GC genomes flag.
- Per tile only for Illumina read-name format; irrelevant for newer patterned flow cells / other platforms.
- Slow on big files (single thread per file, JVM) and memory for long reads; no real streaming of parts of a file.
- No mapping-based information at all: cannot say anything about error rates, alignment, coverage.
- Many panels (11 plots) that most users never open; pass/fail tick on modules that cannot fail (Basic stats).

## 5. Behaviour on long reads

Source: NanoPlot paper (De Coster et al. 2018, Bioinformatics) via search result
https://pmc.ncbi.nlm.nih.gov/articles/pmid/29547981 : Illumina-oriented tools, FastQC among them, are suboptimal for
long reads (85 to 95% accuracy vs <1% error); NanoPlot offers read length histogram (log option), cumulative yield,
length vs quality, and identity/MAPQ bivariates. A direct FastQC-on-ONT benchmark was NOT found (UNVERIFIED).

Reasoned consequences (UNVERIFIED by experiment, follow from the module definitions above):
- Per base quality: a position axis out to tens of kb is meaningless; Q distribution is a property of the read, not
  the position. Median Q ~10-20 for ONT R9 will trip the warn/fail rules (<25 / <20) on perfectly good data
  (modern R10.4.1 Q20+ data passes).
- Per sequence quality and mean Q over 27: ONT/PacBio CLR fail by design; HiFi passes.
- Length distribution: always "warn" since variable; a linear histogram with a 1 Mb tail is unreadable. No N50,
  no log axis, no cumulative yield.
- Duplicates (50 bp prefix) and overrepresented sequences are weak signals; per-base ATGC content and kmer module are not
  informative. Adapter content searches short adapters at read ends only.
- Net: long-read people use NanoPlot/NanoStat/pycoQC/LongQC/SeqKit stats instead; FastQC survives there only as
  habit and via MultiQC.

## 6. Design lessons for SAMStat 2.4

Keep from FastQC:
- A fixed, small plot set, identical order, stable names, one self-contained file (offline, no CDN: this is a
  FastQC strength that SAMStat currently lacks).
- A traffic-light summary at the top with a one-line reason per flag; zero required options.
- Machine-readable summary (JSON + MultiQC already exist) with the same flag IDs as the HTML.
- Short "how to read this" text per plot (collapsible, 2 sentences, what good vs bad looks like).

Avoid or change:
- Do not copy the 12 modules. Candidate standard set (about 6 panels), identical for short and long reads:
  1. Summary table (reads, bases, mean/median length, N50, mean Q, % mapped, %GC) with traffic lights.
  2. Read length histogram, log-x, with N50 marked (this replaces length distribution and gives long-read value).
  3. Per-read mean quality histogram (the FastQC module 3 equivalent) .
  4. Mapping quality / mapped fraction breakdown (SAMStat's existing 4 bins).
  5. Error profile along the read: mismatch, insertion, deletion by position, per MAPQ bin (SAMStat's unique plot);
     for long reads use relative position (0-100%) or fixed bins on a log/binned axis, never raw base positions.
  6. Per-read identity histogram (mapped reads) for long reads - one more plot, but the single most informative
     long-read number; omit when no mapping/MD info.
- Drop: kmer, per-tile, per-base N, overrepresented, adapter and GC-vs-normal modules (low interpretive value,
  platform-specific, or cause false FAILs). If adapter/overrepresentation is wanted, one short table at most.
- Flags: fewer, platform-aware. FastQC's fixed Illumina thresholds mis-flag; in SAMStat pick thresholds only for
  things that are actually wrong regardless of platform (e.g. mapped fraction very low, error rate in the first bin
  extreme, empty file) and show "info" (neutral) instead of fail for the rest. Tolerate that a flag is advisory.
- Subsample for speed, as Falco and `-p` do: a fixed default read cap for plots (documented, reported in the header)
  with exact totals from a cheap counting pass; avoid FastQC's first-100k bias by sampling across the file where possible.
- Never emit per-position series longer than a few hundred points (fixes the 100 MB JSON problem); bin to fixed
  count (e.g. 100 relative bins) for long reads.
- Output compatibility is a distribution asset: keep MultiQC integration and consider matching the fastqc_data.txt
  column names only if it is cheap (Falco's path); otherwise native MultiQC module is the better investment.

## 7. Sources

- FastQC module help: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/
- Falco: https://github.com/smithlabcode/falco
- NanoPlot / NanoPack paper: https://pmc.ncbi.nlm.nih.gov/articles/pmid/29547981
