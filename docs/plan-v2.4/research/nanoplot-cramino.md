# Long-read QC tools: what matters (NanoPlot, cramino, chopper, nanoq, pycoQC, ToulligQC, LongQC)

Researched 2026-10-14. Items without a URL or marked UNVERIFIED are from general knowledge or inference, not checked.

## 1. Tool by tool

| Tool | Input | Metrics / plots | Speed / status | Source |
|---|---|---|---|---|
| NanoPlot (v1.48.0 on the releases page) | FASTQ, FASTA, sequencing_summary, sorted BAM, uBAM, CRAM, feather/arrow | Read-length histogram with N50 marked; log-transformed length histogram (`--loglength`); length vs basecall quality bivariate; yield over time; channel heatmap; length/quality violins over time. BAM only: aligned vs sequenced length, identity vs length, identity vs quality, MAPQ vs length. NanoStats text/TSV table. | Summary-file input is fastest; kde plots slow; downsampling happens after collection so saves little. Python. Optional `--include-js embedded` makes the report offline but grows a full run from <1 MB to ~75 MB. | https://github.com/wdecoster/NanoPlot , https://github.com/wdecoster/NanoPlot/releases |
| cramino | BAM/CRAM, stdin, remote (http/ftp/s3) | Text/JSON/TSV table only: reads, yield, yield >25 kb, N50, N75, median/mean length, aligned N50/N75, median/mean/modal identity (gap-compressed, from `de` or NM+CIGAR; "estimated identity" from base quality for uBAM). Optional `--hist` length/identity histograms, `--arrow` feeding NanoPlot/NanoComp. | Claims 140 Gbase BAM in 12 min, <1 GB memory. Rust. | https://github.com/wdecoster/cramino |
| chopper | FASTQ | Filter/trim only (mean Q, min/max length, GC, contaminants). No report. | Rust rewrite of NanoFilt; no quantified claim; runtime often dominated by gzip decompression. | https://github.com/wdecoster/chopper |
| nanoq | FASTQ/FASTA | Summary stats and filtering. Uses nanopore-appropriate mean quality (probability-space, not mean of Phred). | Own docs claim up to ~450x vs NanoStat in fast mode (no quality); ~1.3-1.5x vs seqtk fqchk on large data; 5-10x less memory. Rust. | https://joss.theoj.org/papers/10.21105/joss.02991.pdf , https://docs.rs/nanoq/ |
| pycoQC v2.5.2 | sequencing_summary.txt (Albacore/Guppy) | Length and quality distributions, length vs quality, channel activity, yield/quality/length over time, alignment plots. Interactive HTML + JSON. | Archived 2024-10-18, unmaintained. | https://github.com/a-slide/pycoQC |
| ToulligQC | Guppy/Dorado sequencing summary + telemetry | Plots, text stats, HTML report; barcode support; Nextflow and Snakemake wrappers. | Active; run-level, ONT only. | https://github.com/GenomiqueENS/toulligqc |
| LongQC | FASTA/FASTQ/PacBio subread BAM | Length, GC, coverage dispersion, end-adapter patterns, sdust low complexity, reference-free via self-mapping; platform QC for PacBio/ONT. | Heavier (minimap2 self-mapping, Docker recommended); 31 open issues. | https://github.com/yfukasawa/LongQC |

Could not fetch (404 or no page): nanoq and ToulligQC GitHub READMEs at the guessed URLs; used search summaries instead. Exact nanoq/NanoPlot version dates UNVERIFIED (releases page showed no year).

## 2. Observations

1. Two families. Run-level tools (pycoQC, ToulligQC) need ONT sequencing_summary and are about the flowcell (channels, time). They do not apply to SAM/BAM/FASTQ files from arbitrary sources. SAMStat should not go there.
2. File-level tools (NanoPlot, cramino, nanoq) answer the same few questions: how much, how long (N50), how good (Q or identity). cramino shows a text table of ~12 numbers is enough for many users and is the fastest path.
3. The headline numbers are a fixed list: reads, total bases, N50, mean/median length, mean/median Q (or identity), and yield above thresholds. NanoPlot's two core plots are the length histogram (linear and log) and length vs quality.
4. Log-scale length is essential for long reads: length distributions are heavy-tailed (ultra-long libraries span 1 kb to >1 Mb), so a linear histogram hides the bulk. Mark N50 on it.
5. Identity (BAM only) is the long-read counterpart of SAMStat's mismatch/indel profiles, and is the metric cramino and NanoPlot lead with. Gap-compressed identity vs blast-style identity must be defined explicitly in the report (cramino documents its choice).
6. Quality: the correct mean read Q averages error probabilities, not Phred values (nanoq). Averaging Phred overstates quality. SAMStat's cap at Q41 is wrong for modern long reads (Q50+ for HiFi; UNVERIFIED for exact ranges) and should be revisited.
7. Speed expectations are set by Rust tools: cramino ~140 Gbase/12 min; nanoq summaries are bound by decompression. Python NanoPlot is slow on large BAMs (UNVERIFIED benchmark; qualitative from its own docs recommending summary files and subsampling).
8. Offline reports: NanoPlot itself added embedded plotly because CDN reports fail offline, and pays ~4.8 MB per HTML file. This confirms the SAMStat CDN weakness, and suggests static SVG or tiny inline JS instead of embedding Plotly.
9. What users like (UNVERIFIED, general experience and inferred from design): NanoPlot's standard plot set and one-line CLI; cramino's speed and single table; nanoq's speed. Dislikes (UNVERIFIED, plausible from the docs): NanoPlot slow with large inputs and many output files; kde plots slow; hex plots break with the plotly path; pycoQC abandoned; LongQC heavy install.

## 3. Minimum set that matters (recommendation for SAMStat 2.4)

Summary table (always, at top of report and in JSON/MultiQC): reads, bases, mean/median length, N50, (N90 optional), mean Q (probability-space), % reads Q>=10/Q>=20 (or equivalent), and for aligned data: % mapped, median identity (gap-compressed), per-MAPQ-bin counts.

Plots (fixed set, same every time):
1. Read length histogram, log-x bins, N50 marked (replaces per-position plots as the primary long-read view).
2. Mean read quality histogram.
3. Identity histogram (aligned input only).
4. Existing SAMStat per-position mismatch / insertion / deletion / quality plots, but binned to a fixed number of positions (e.g. <=100 bins on a relative or log axis) so a 1 Mb read cannot produce 100 MB output.
5. Keep MAPQ stratification (SAMStat's differentiator; NanoPlot only offers MAPQ scatter).

Do not add: time/channel/flowcell plots (need sequencing_summary, run-level), bivariate kde/hex scatter plots (slow, hard to read at a glance; a 2D length-vs-quality heatmap could be considered later but is not minimum), GC/low-complexity, contaminant filtering (chopper's job), trimming/filtering.

Performance targets to state: stream input, no per-read storage (N50 needs only a length histogram with fine log bins, or a sorted length array of 8 bytes/read), stdin support (cramino does it), parallel htslib decompression threads for a single BAM (`hts_set_threads`) to approach cramino's throughput.

Presentation: auto-detect long vs short from length distribution (e.g. median > 1 kb or max > ~1 kb) and switch to log length axis; no required flag.

## 4. Gaps this leaves open (for the planner)

- Exact definition of identity to use (gap-compressed like cramino/minimap2 `de`, or matches/alignment-columns).
- Whether to read `de`/NM tags instead of MD (MD absent in many long-read BAMs; UNVERIFIED how SAMStat behaves when MD is missing, check src/sambamparse).
- Whether 2D length-vs-quality heatmap earns a place.
