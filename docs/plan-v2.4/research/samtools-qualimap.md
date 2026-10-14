# Alignment-QC landscape: samtools stats, Qualimap, Picard, mosdepth, alfred, BamQC

Date: 2026-10-14. Verified = fetched this session. Everything else is marked UNVERIFIED (from memory).

## 1. Tool summary

| Tool | Metrics / plots | Speed / threading | Main use |
|---|---|---|---|
| samtools stats + plot-bamstats | Verified sections: SN summary numbers, FFQ/LFQ per-cycle quality, GCF/GCL, GCC/FBC base content per cycle, IS insert size, RL read length, MAPQ distribution, ID indel size, IC indels per cycle, COV, GCD, MPC mismatches per cycle x quality (reference required). plot-bamstats makes static PNGs through gnuplot (UNVERIFIED). | `-@` adds compression threads only (default 0), so it is effectively single-core (verified). C, streaming, fast. | De facto cheap text summary. MultiQC reads it. Plots are dated. |
| Qualimap bamqc | UNVERIFIED (doc page failed to load, SSL error): coverage across reference, coverage histogram, GC vs coverage, mapping quality across reference and histogram, insert size, clipping profile, mismatch/indel summary, homopolymer indels. HTML + PDF. | Java, `-nt` threads, `--java-mem-size`. Reputed slow and memory-hungry on big BAMs (UNVERIFIED). | WGS/exome coverage uniformity. Has a long-read-unfriendly design (UNVERIFIED). |
| Picard CollectAlignmentSummaryMetrics | Verified columns: PF_MISMATCH_RATE, PF_HQ_ERROR_RATE, PF_INDEL_RATE, STRAND_BALANCE, PCT_CHIMERAS, PCT_ADAPTER, MEAN_READ_LENGTH, PF_HQ_ALIGNED_Q20_BASES, BAD_CYCLES, in-pair counts. Per category (FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR). Mismatch rates need a reference (inferred). | Java, single-threaded (UNVERIFIED). Text table only; plots via separate collectors (insert size histogram, quality by cycle). | Pipeline gate metrics (GATK ecosystem). Numbers, not figures. |
| mosdepth | Verified: depth only. Per-base, per-window or per-region mean, cumulative distribution, threshold counts, quantized bins. No built-in plotting beyond `plot-dist.py` giving `dist.html`. | Verified: Nim, ~2x faster than `samtools depth` (30X genome 25 min single-threaded per-base vs 50 min); `-t` decompression threads, flattens beyond ~4. | Coverage, fast. Complements, does not compete with, error QC. |
| alfred (qc) | Verified only: "interactive multi-sample BAM alignment statistics", short and long reads. Details at https://www.gear-genomics.com/docs/alfred/ not read. UNVERIFIED from memory: per-read-group read length, MAPQ, base quality, insert size, GC, coverage, error rates, JSON output and a web viewer. | C++/htslib, UNVERIFIED threading. | The closest feature overlap with a "short + long" aligned QC tool. |
| BamQC (Babraham) | UNVERIFIED: FastQC-like for BAM, from Babraham. I could not confirm its current status or features; may not be maintained. Treat as not a live competitor. | UNVERIFIED | - |

Sources: https://www.htslib.org/doc/samtools-stats.html ; https://github.com/brentp/mosdepth ; https://broadinstitute.github.io/picard/picard-metric-definitions.html ; https://github.com/tobiasrausch/alfred

## 2. What nobody covers well

1. **MAPQ-stratified error profiles.** samtools MPC gives mismatches per cycle x base quality, not per MAPQ class. Picard splits only by read category and reports one scalar rate (PF_HQ_ERROR_RATE uses a MAPQ>=20 high-quality subset, UNVERIFIED detail). Qualimap and mosdepth have no per-position error profile. Nothing found shows mismatch/insertion/deletion along the read for MAPQ strata side by side. This is SAMStat's niche.
2. **Reference-free for the user.** samtools MPC and Picard need a reference FASTA. SAMStat uses the MD tag and needs none (limitation: MD must be present).
3. **FASTA/FASTQ and SAM/BAM/CRAM in one tool with one report.** The others are alignment-only. FastQC covers FASTQ and not alignments.
4. **One self-contained standard page.** samtools yields text plus gnuplot PNGs, Picard text, Qualimap a heavyweight Java report, mosdepth a distribution only. No tool offers "the same fixed plots every time" for aligned reads, FastQC style.
5. **Long reads.** samtools and Picard produce tables that stay usable. Per-cycle plots are not meaningful past a few kb. Alfred claims long-read support (verified claim only). Gap: read-length histogram / N50 / identity are not part of the classic suite (SAMStat shares this gap).

## 3. Where SAMStat fits

- Position it as the "FastQC for the mapped file": run after the aligner, answers "is the alignment good, and where do errors come from, and do low-MAPQ reads look different?"
- Do not chase coverage (mosdepth), insert-size/GC-bias depth models or pipeline gate tables (Picard). Keep coverage out.
- Speed benchmark to beat: samtools stats on one BAM, with similar single-core throughput. SAMStat's one-file-per-thread model gives zero gain for a single big BAM. A reader-thread plus worker-thread scheme (or htslib `-@` decompression threads at minimum) is the cheapest fix. The competitors' own threading is also modest (samtools `-@` compression only; mosdepth plateaus at ~4), so decompression threads plus one parse thread would reach parity or better.

## 4. Recommendations (simple, FastQC-like)

1. Keep the fixed plot set: summary table, MAPQ distribution, read length, base quality along read, mismatch / insertion / deletion along read, each stratified by MAPQ class. Same plots every run.
2. Add for long reads only the minimum: log-scale read-length histogram, N50 in the summary, and binned (not per-position) profiles past ~1 kb. No new plot families.
3. Emit samtools-stats-compatible numbers in the summary (names of SN fields) where cheap, so users can cross-check. Do not duplicate the SN text format.
4. Ensure the MultiQC output carries the MAPQ-stratified profiles, since that is the unique content others lack.
5. Add htslib decompression threads (`hts_set_threads`) before any larger threading rewrite.
6. Skip: coverage, GC-bias, insert-size modelling, PDF output, a GUI.

## 5. Unverified items to check before publishing claims
Qualimap feature list and speed (https://qualimap.conesalab.org/ failed to load); alfred qc output list; BamQC status; Picard threading; plot-bamstats dependence on gnuplot.
