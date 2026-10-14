# Competitor speed survey and SAMStat throughput target

Date: 2026-10-14. Method: web search only; primary papers on PMC were blocked by a captcha on fetch, so several numbers come from search-result snippets and are flagged. No local benchmarking was done.

## 1. Published numbers

| Tool | Input | Figure | Source | Status |
|---|---|---|---|---|
| cramino | 140 Gbase long-read BAM | about 12 min, default 4 threads, under 1 GB RAM (= about 11.7 Gbase/min, about 195 Mbase/s) | https://biostars.org/p/9544047 , https://pmc.ncbi.nlm.nih.gov/articles/PMC10196664 | Vendor claim in tool docs/announcement; not independently verified |
| cramino vs NanoPlot/NanoStat | BAM | Paper positions cramino as replacement for NanoStat because it is "too slow" for current yields. No head-to-head table found | PMC10196664 | UNVERIFIED (numbers) |
| Falco vs FastQC | uncompressed FASTQ | Test 1: 48 s vs 3 m 30 s; test 5: 15 m 49 s vs 41 m 27 s (about 2.6-4.4x) | https://pmc.ncbi.nlm.nih.gov/articles/PMC7845152/table/T3 | From search snippet; file sizes not captured |
| Falco vs FastQC | gzip FASTQ | Test 1: 1 m 19 s vs 3 m 49 s; test 9: 2 s vs 7 s (about 3x) | https://pmc.ncbi.nlm.nih.gov/articles/PMC7845152/table/T4 | Same caveat |
| fastp | FASTQ | Same order as FastQC, slower than Falco on most files in the Falco paper (fastp does more: trimming) | PMC7845152 tables 3-4 | Qualitative only; no fastp absolute numbers verified |
| mosdepth | about 30x human WGS BAM | 25 m 23 s at 1 thread, 14 m 27 s at 3 threads, about 1.2 GB RAM; no gain beyond about 4 decompression threads. README claims <9 min in window mode at 3 threads | https://pmc.ncbi.nlm.nih.gov/articles/PMC6030888 | Paper (2018); README figure is secondary |
| samtools depth | same 30x BAM | 50 m 12 s, 1 thread (from mosdepth paper table) | PMC6030888 | Not `samtools stats`; proxy only |
| samtools stats | WGS BAM | No published figure found. Single-threaded parse of every record, so expected near samtools depth class (tens of minutes at 30x) | - | UNVERIFIED |
| Qualimap bamqc | 30x BAM | No figure found. JVM, `-nt` threads, needs large `--java-mem-size`. Widely reported as slow (tens of minutes to hours) | - | UNVERIFIED |
| NanoPlot | 100 Gb ONT | No figure found. Python/pandas, plots every read; commonly reported hours for large runs | - | UNVERIFIED |
| FastQC | FASTQ | Single-threaded per file, JVM; about 3-4x slower than Falco per above | PMC7845152 | Derived |

## 2. What "fast" means in 2026

The reference points are:
- Short-read FASTQ QC: Falco/FastQC finish a typical file in minutes; users tolerate this because it runs once per sample.
- Alignment QC: mosdepth at about 15 min for 30x is the "fast" bar; samtools-class tools take 30-50 min.
- Long-read QC: cramino at about 12 min per 140 Gbase on 4 threads is the bar that displaced NanoStat/NanoPlot for yield-scale data.

Implied throughput (my arithmetic from the figures above, assuming 30x human = about 90 Gbase):
- cramino: about 195 Mbase/s on 4 threads.
- mosdepth 3 threads: about 105 Mbase/s (90 Gbase / 867 s).
- Falco and FastQC are bound by gzip decompression (about 100-300 MB/s per stream) and are in the same range or lower.

## 3. Proposed SAMStat targets (recommendation, not a verified fact)

1. Primary bar: 30x human short-read BAM (about 90 Gbase) in at most 10 min on 8 threads, and no slower than 15 min on 4 threads (match mosdepth).
2. Long reads: 100 Gbase ONT BAM in at most 10 min on 4-8 threads (match or beat cramino's 12 min / 140 Gbase on 4 threads, i.e. at least 190 Mbase/s).
3. FASTQ/FASTQ.gz: at least Falco parity (about 3x FastQC); gzip bound, so a parallel decompressor (e.g. libdeflate/bgzf threads) matters more than the stats code.
4. Memory under 2 GB regardless of read length or file size (cramino under 1 GB; mosdepth 1.2 GB).
5. Output size independent of read length (a 1 Mb read must not produce 100 MB JSON): fixed-size binned summaries.

## 4. Implications for SAMStat 2.4

- The decisive gap is single-core-per-file. htslib decompression threads (hts_set_threads, `-@`) are the cheapest win; the mosdepth data show returns diminish beyond about 4 decompression threads, so real scaling needs parallel record processing (per-thread accumulators merged at the end, which suits SAMStat's additive histograms) or region/BGZF-block splitting for indexed files.
- CRAM and BAM decode cost dominates; the MD/CIGAR walk adds per-base work that cramino and mosdepth avoid, so bound per-base work on long reads (bin positions, cap per-position resolution at a few kb).
- Do not compete with NanoPlot on plot count; compete on wall time and on one fixed set of plots.
- Measure, don't assume: add a reproducible benchmark script (fixed 30x short BAM, fixed ONT BAM) to CI-adjacent docs so the targets above are checked on SAMStat itself.

## 5. Gaps

Not verified: NanoPlot, Qualimap, samtools stats, fastp absolute timings; hardware used in the Falco paper; cramino independent benchmark. Re-check Falco/fastp tables directly at https://pmc.ncbi.nlm.nih.gov/articles/PMC7845152 and cramino at https://github.com/wdecoster/cramino.
