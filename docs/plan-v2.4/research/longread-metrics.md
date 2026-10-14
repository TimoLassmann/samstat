# Long-read alignment QC metrics for SAMStat 2.4 (ONT R10.4.1 / Dorado, PacBio HiFi / Revio)

Date: 2026-10-14. Scope: which metrics tell a good run from a bad one at a glance, and the cheapest way to compute them in SAMStat's one-pass parser. Items without a URL come from general knowledge and are marked UNVERIFIED where relevant.

## 1. Verified external facts

| Fact | Source |
|---|---|
| minimap2 PAF col10/col11 = BLAST-like identity; `NM` = mismatches + gaps (counted per gap base); `de` = gap-compressed per-base divergence; `dv` = approximate divergence (no formulas given on that page) | https://lh3.github.io/minimap2/minimap2.html |
| Dorado aligns with minimap2, default preset `lr:hq`, overridable by `--mm2-opts` | https://github.com/nanoporetech/dorado/blob/release-v1.0/README.md |
| Dorado documents tags `dx:i` (duplex 1 / 0 / -1), `TS:A`, `RX:Z`, `BC`, `RG`, `pt:i` (poly-A), `mv` (only with `--emit-moves`). README does NOT document `qs`; `qs` appears only in an example FASTQ header. | same |
| PacBio BAM: `rq` float in [0,1] = predicted read accuracy; `np` = number of complete passes; QUAL is Phred. README does not describe QUAL binning. | https://pacbiofileformats.readthedocs.io/en/latest/BAM.html |
| NanoPlot plot set: read length hist (+N50), log-length hist, length vs quality, yield over time, channel heatmap, and (BAM) aligned vs sequenced length, identity vs length, identity vs quality, MAPQ vs length | https://github.com/wdecoster/NanoPlot |

UNVERIFIED (from memory, confirm before implementing): ONT `qs` tag = mean read Q-score (float) written by Dorado BAM output, and `rq` is also written by pbmm2/ccs as the HiFi predicted accuracy. PacBio Revio/Sequel IIe HiFi base QUAL is binned to a small set of values (reported as Q0/Q10/Q20/Q30/Q40-ish bins, i.e. not a 0-93 continuum); check actual BAMs. Dorado simplex SUP R10.4.1 median read accuracy about 99% (Q20) and duplex about Q30: vendor claims, UNVERIFIED.

## 2. What distinguishes good from bad (ranked by signal per plot)

1. **Read-length distribution, log10 x-axis, with N50 and mean/median stated.** One plot answers: right fragment size? degraded/short? ultra-long tail? Linear axis is useless for ONT (heavy tail). HiFi: expect a tight peak (typically 10-25 kb). Bad run: bimodal with a short-fragment peak (adapter dimers, degraded DNA), or a peak far from the library target.
2. **Per-read identity histogram (0.5-1.0, or Q-scaled), one curve per MAPQ class or just primary.** The single most informative number for chemistry/basecaller health. Good: ONT SUP R10.4.1 mode around 0.98-0.99 (UNVERIFIED), HiFi mode >= 0.999 (UNVERIFIED but consistent with rq >= 0.99 definition of HiFi, https://pacbiofileformats.readthedocs.io/en/latest/BAM.html semantics). Bad: wide left shoulder or second mode below 0.9 (contamination, wrong reference, failed pore chemistry).
3. **Error-type breakdown: mismatch / insertion / deletion rates (per aligned base), plus homopolymer-indel fraction.** Distinguishes platform issues: older ONT is deletion/homopolymer heavy, HiFi residual errors are mainly indels in long homopolymers, basecaller-model mismatch shows as elevated mismatches. Plain three-bar rate is enough; homopolymer split is optional (see Section 4).
4. **Mapped fraction and MAPQ class (existing 4 bins) plus supplementary and secondary fractions.** Bad: mapped < ~90% on a known genome (contamination/wrong reference); supplementary rate high suggests chimeras or SVs/library artefacts (library prep ligation chimeras are an ONT issue, UNVERIFIED rate).
5. **Aligned fraction of read** (aligned query bases / read length, soft+hard clip excluded). Catches adapters/chimeras/junk ends. Cheap, already derivable from CIGAR.
6. **Read Q-score distribution (mean per read) vs empirical identity.** Calibration check: the mean per-read Q should roughly agree with identity-derived Q. Useful for ONT `qs`, HiFi `rq`. Lower priority than 1-3; fold into one sentence of summary or a scatter only if it is cheap.

Not worth adding: yield-over-time, channel heatmap (need run metadata not in BAM, ONT-specific, NanoPlot/MinKNOW cover them), kinetics tags, per-position plots beyond a few kb.

## 3. Definitions (exactly one each; document in the report)

- **Per-alignment identity (BLAST-like)** = (M_bases - mismatches) / (M + I + D) where M counts CIGAR M/=/X, I and D count inserted/deleted bases. Equals minimap2 col10/col11.
- **Gap-compressed identity** = 1 - (mismatches + number of indel *events*) / (M + number of indel events) (equivalent to 1 - `de` in minimap2 terms, https://lh3.github.io/minimap2/minimap2.html). Better than BLAST identity for long reads because one 30 bp deletion does not count as 30 errors; use this as the headline "read identity". Compute from CIGAR plus NM (mismatches = NM - I - D bases) so the MD tag is not required. If `=`/`X` CIGAR ops are present use them directly.
- Report identity per **primary** alignment, weighted by one alignment per read; for chimeric reads add the supplementary only to the supplementary count (avoid double counting reads in the histogram). Alternatively length-weight the identity (base-weighted) in the summary line; histogram stays per-read.
- **N50** = length L such that reads >= L contain half of all bases. Compute from the length histogram on read length (query length incl. hard-clip reconstructed from CIGAR for supplementary rows is not needed: use primary only).
- **Homopolymer indel**: indel adjacent to a run of >= 4 identical reference bases (threshold UNVERIFIED as convention; pick one and print it). Requires the reference or MD; without reference only the read-side run is available. This is the one metric that needs reference access, so default OFF/omitted unless the reference is given. Recommendation: do not implement in 2.4 core.

## 4. Fit with the current code (2.3.0) and cost

Existing: `parse_alignment` already derives CIGAR + MD per-position mismatch/ins/del; per-position plots are interval-binned. Gaps for long reads:

| Need | Where | Cost | Note |
|---|---|---|---|
| Read length histogram (log bins, e.g. 10 bins/decade from 100 bp to 10 Mb) + N50 | collect.c (new fixed-size array, no per-position storage) | small | Fixed memory: ~60 counters. Replaces unbounded per-position growth. |
| Per-read identity histogram | collect.c: accumulate M, I, D events, mismatches from CIGAR/NM; bin 0.50-1.00 in 0.5% steps, or finer (0.1%) in 0.95-1.00 for HiFi | small | Use NM when present, else MD. Two scales via one array of ~1000 bins from 0 to 1. |
| Error-rate totals per MAPQ class | collect.c | trivial | Already half there (counts per position). |
| Supplementary/secondary fractions | flag test | trivial | 2.3 tests already cover counts. |
| Aligned fraction | CIGAR clip ops | trivial | |
| Per-position plots for long reads | plot | cap at first N bases (e.g. 5 kb?) or switch to relative-position (percent of read) x-axis | Root cause of the ~100 MB JSON for a 1 Mb read: unbounded position arrays. Fix at the collector (bin by position, cap), not the JSON writer. |
| Q-score tags (`qs`, `rq`) | optional tag reads | trivial | If tag absent, derive mean Q from QUAL (needs mean of error probabilities, not mean of Q). |

Key point: the genuinely new plots (length, identity) are fixed-size histograms, so they add almost no memory and make single-file runtime unaffected. Per-read identity needs NM or MD; minimap2 with `--MD`/default NM usually supplies NM; pbmm2 supplies NM (UNVERIFIED). If neither tag exists, identity panel prints "n/a (no NM/MD)" rather than failing.

## 5. Recommended "standard set" for long reads (same page layout for short reads where meaningful)

Keep a FastQC-like fixed order, identical every run:

1. Summary table: reads, bases, mean/median length, **N50**, mapped %, supplementary %, secondary %, median identity (gap-compressed), median read Q.
2. Read length histogram, log x, N50 marker.
3. Per-read identity histogram (primary, split by MAPQ class using the existing bins; keep colours stable).
4. Error profile: mismatch / insertion / deletion rate per MAPQ class (existing SAMStat strength), per position for short reads, per relative position (or first-N-bases) for long reads.
5. Base-quality (existing) plus mean-Q histogram.

Auto-detect mode with no option: if median read length > ~1 kb (or max > ~10 kb), long-read layout (log-length and identity first, per-position capped); otherwise short-read layout unchanged. Include a one-line "reading guide" caption under each plot with healthy ranges (example: "HiFi: peak >= 99.9%"). Ranges are guidance, label as typical, not pass/fail, to avoid maintenance churn as chemistry evolves.

## 6. Decisions requiring judgement (for the planning lead)

- Headline identity: gap-compressed (recommended) vs BLAST-like. Show only one in plots; list both in the JSON.
- Do NOT add pass/fail thresholds in 2.4; show captions with typical values only.
- Defer homopolymer analysis, kinetics, time-series, channel maps.
- Add JSON/MultiQC fields: n50, median identity, identity histogram, length histogram (fixed-size arrays, so schema stays small).

## 7. Open checks before implementation

- Confirm real-world tags in a Dorado R10.4.1 BAM (`qs`, `du`, `ns`, `NM`, `MD` presence) and a Revio pbmm2 BAM (`rq`, `np`, QUAL binning). Obtain test BAMs and measure.
- Confirm cost of identity computation on a 1 Mb read (should be O(CIGAR ops)).
