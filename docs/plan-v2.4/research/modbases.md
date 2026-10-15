# Modified-base (MM/ML) summaries for SAMStat 2.4

Researched 2026-10-14. Source tag: [V] = verified from the cited URL, UNVERIFIED = not confirmed.

## 1. Recommendation (short)

- Do NOT add a modbase section to the default report. FastQC-style "same plots every time" is better served by keeping the standard set fixed.
- DO add a tiny, conditional **"Base modifications" block** that appears only if at least one record carries MM/ML. No option, no plots beyond one. If no MM/ML is present, output is byte-identical to today (goldens unchanged).
- Minimal useful version = 5 numbers/tables, all computed from one pass over reads that SAMStat already parses:
  1. % of reads with an MM tag (and with a non-empty MM) out of all primary reads.
  2. Per canonical base + mod code (e.g. C/m, C/h, A/a): number of calls (all, and at p >= 0.5).
  3. Fraction called modified at p >= 0.5 per code (one fixed threshold, not tunable).
  4. One probability histogram (ML value, 16 or 32 bins of the 0-255 range) per mod code. This shows whether the model is confident (bimodal) or not.
  5. Count of reads with MM/ML inconsistent with SEQ (MN mismatch, or coordinates past sequence end) - catches hard-clipped supplementary alignments (see section 5).
- Optional single plot: mod fraction (at p >= 0.5) by MAPQ bin, reusing SAMStat's existing four MAPQ bins. This is the one place SAMStat is distinctive; but defer it until users ask. UNVERIFIED demand.

Verdict on value: moderate and growing for ONT/PacBio, because modBAMs are now the default output of both platforms and "did my BAM keep its mod tags?" is a real failure mode (hard clipping, tools dropping tags). A presence/consistency check is high value, cheap, and fits QC. Rich analysis (pileups, motifs, per-site) is out of scope - that is modkit / pb-CpG-tools territory.

## 2. What other tools report

### modkit (Oxford Nanopore)
- Repo: https://github.com/nanoporetech/modkit ; latest release v0.6.4, published 2026-06-11 [V] (https://api.github.com/repos/nanoporetech/modkit/releases/latest).
- `modkit pileup`: modBAM -> bedMethyl with 18 columns incl. Nvalid_cov, fraction modified, Nmod, Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall [V]. Per-site, needs sorted+indexed BAM. Not QC.
- `modkit modbam summary` [V] (https://github.com/nanoporetech/modkit/blob/master/book/src/intro_summary.md; the older `modkit summary` name appears in the README):
  - Header (`#` lines): bases with calls, total_reads_used, count_reads_{base}, modification codes per base, filter threshold per base.
  - Table columns: base, code (`-` = canonical), pass_count, pass_frac, all_count, all_frac. "pass" = call probability above a threshold, chosen automatically as the 10th percentile of calls (`--filter-percentile`; default not stated on the page, example log shows 10) or fixed via `--filter-threshold`.
  - Uses the whole BAM by default since v0.6.0; `--num-reads N` / `--sampling-frac F` to subsample; `--region`, `--motif`, `--matched-only` exist.
  - Fractions are fractions of CALLS, not of reads. No fraction-of-reads-with-MM headline is printed (only reads used / per base).
  - Probability histograms: `modkit modbam sample-probs` produces them (name from the docs page above; histogram output details UNVERIFIED).
- Takeaway: modkit's QC = per-code call counts + pass fraction + a threshold. That is exactly the minimal set above. Its adaptive percentile threshold is a complexity we should skip.

### pb-CpG-tools (PacBio)
- Repo: https://github.com/PacificBiosciences/pb-CpG-tools [V].
- Input: aligned, indexed BAM/CRAM with MM/ML; only the `C+m` key is processed; other keys (e.g. `A+a`) and combined keys like `C+mh` are ignored/unsupported [V].
- Output: combined (and hap1/hap2 if HP tag) bed.gz + bigWig of per-CpG modification score and coverage [V]. No separate QC/summary report [V].
- Takeaway: downstream site-level tool, not QC. It shows the convention to watch: PacBio uses `C+m` on CpG only.

### Others (not verified)
- ONT wf-basecalling / wf-human-variation report modification summaries via modkit - UNVERIFIED details.
- MultiQC has no modkit module known to me - UNVERIFIED; check before promising integration.
- samtools/ pysam can read tags but do not summarise.

## 3. Format facts (SAM spec)

Source: https://raw.githubusercontent.com/samtools/hts-specs/master/SAMtags.tex [V]
- MM:Z regex `([ACGTUN][-+]([A-Za-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*`. First char = unmodified base, `+/-` = strand relative to as-sequenced SEQ, codes are letters (m, h, a, ...) or ChEBI numbers. Positions are deltas counting only that base type from the 5' end of SEQ. Empty list (`C+m;`) is legal.
- Flag after code: `?` = skipped bases have unknown status; `.` or none = skipped bases assumed low probability (implicit unmodified). Matters: with `?`, "not listed" is NOT evidence of canonical.
- ML:B:C = one byte per listed position, in MM order. Value N means probability in [N/256, (N+1)/256). Absolute probabilities.
- MN:i = SEQ length when MM/ML were written; if it differs from current SEQ length, MM/ML are out of date. MN was added to the spec in Sept 2024 [V]; MM/ML became standard Feb 2022 [V]. Many older BAMs lack MN.
- Reverse-strand (0x10): tag refers to as-sequenced orientation, i.e. opposite to SEQ in the BAM. htslib handles this internally.

## 4. htslib API (hts_base_mod)

Declarations from htslib 1.16 sam.h (https://raw.githubusercontent.com/samtools/htslib/1.16/htslib/sam.h) [V]:
- `hts_base_mod { int modified_base; int canonical_base; int strand; int qual; }` - modified_base = char code or negative ChEBI; canonical_base in A/C/G/T/N; strand 0/1; qual = 256*prob, -1 if unknown.
- `hts_base_mod_state *hts_base_mod_state_alloc()/..._free()`.
- `int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state)` - parse MM/ML, 0 ok / -1 fail.
- `int bam_mods_at_next_pos(b, state, mods, n_mods)` - iterate every query position; returns count (may exceed n_mods).
- `int bam_next_basemod(b, state, mods, n_mods, int *pos)` - skip to next position that HAS mods; 0 = done. Best for a summary: loop cost proportional to number of calls, not read length.
- `int bam_mods_at_qpos(...)` - random-ish access with increasing qpos.
- `int *bam_mods_recorded(state, &ntype)` and `bam_mods_query_type(state, code, &strand, &implicit, &canonical)` - list codes present and the implicit flag.
- The 1.16 header I fetched did not show HTS_MOD_* flag constants [V as observed; the fetch tool may have trimmed them] - UNVERIFIED whether 1.16 has them.

Later changes (https://github.com/samtools/htslib/blob/develop/NEWS) [V]:
- 1.18: MN tag support; API improved so callers can distinguish unchecked from unmodified bases; `bam_mods_queryi()` added.
- 1.19: docs. 1.20: better validity checking of mod tags.
- 1.23: fix segfault on empty valid MM tag; fix `bam_next_basemod` with HTS_MOD_REPORT_UNCHECKED.
- Implication: SAMStat bundles 1.16. Empty-MM segfault (fixed 1.23) and weaker validity checks apply to 1.16. UNVERIFIED whether the empty-MM crash exists in 1.16 specifically. If modbases are adopted, bump bundled htslib (to a current release, >= 1.23) - this is a prerequisite and has its own regression risk (goldens, packaging tests, bioconda recipe). Check how the plan treats an htslib bump before committing.
- Repo check: `grep` shows SAMStat currently has no hts_base_mod usage.

## 5. Pitfalls specific to SAMStat

1. Supplementary/hard-clipped alignments: MM/ML positions are relative to the original read; after hard clipping SEQ is shorter, so tags are invalid unless MN-aware. Count and exclude those from mod stats; report the count (item 5 in section 1). Without MN, fall back to checking positions against SEQ length.
2. Primary only: count mods on primary alignments (not secondary/supplementary) to avoid double counting. SAMStat already tracks secondary/supplementary counts (2.3 tests).
3. Unmapped BAM (uBAM, direct from basecaller) is the common input for mod QC; SAMStat handles unmapped reads, so this works without reference. Do not require alignment (unlike pileup tools).
4. CRAM: tags survive; no special handling.
5. Cost: parsing MM/ML for every read adds time. ONT human 5mC reads have ~1 call per CpG (about 1 per 100 bp); `bam_next_basemod` keeps this cheap, but long reads at Mb scale produce 10^4-10^5 calls per read. Benchmark before accepting; -p sampling already bounds it. Measured numbers: none (not run).
6. Memory/output: store only fixed-size histograms per (base, code, strand), never per-read or per-position data. This avoids the JSON-size problem seen with 1 Mb reads. JSON: add one small `modifications` object; MultiQC: one small table (+ one bar/histogram). Schema bump needed (additive).
7. Mod codes are open-ended (ChEBI numbers, `N` base). Cap distinct (base,code) keys (e.g. 16) and lump the rest as "other" to keep the output bounded.
8. Strand: report mods per canonical base/code, merge strands in the headline; keep strand in JSON only if free.

## 6. Minimal design sketch (if adopted)

Per read (primary, MM present): `bam_parse_basemod`; loop `bam_next_basemod` with an array of 16 `hts_base_mod`; for each call update `hist[key][qual>>4]` (16 bins), `calls[key]++`, `mod[key] += (qual >= 128)`. Count reads with MM, reads with empty MM, reads failed/inconsistent.
Report block (HTML): one small table (code, calls, % modified at p>=0.5) and one histogram plot with one line per code. Shown only when reads_with_MM > 0. No new CLI flags.
Effort estimate: roughly 200-300 lines C plus tests and fixtures (need a small modBAM fixture; make with samtools/pysam script, UNVERIFIED availability) plus schema/doc updates. Highest cost is the htslib bump, not the feature.

## 7. Priority call for the maintainer

- Rank: below long-read length/N50/identity work (core to "long read QC"), above any new input formats.
- Cut line: if schedule is tight, ship only items 1 and 5 (presence and consistency of MM/ML as two counts). That already prevents the worst real-world failure (silently lost or invalid mod tags) at nearly zero cost and zero new plots.
- Do not implement: adaptive thresholds, motif/region filters, per-site pileups, haplotype splits, per-position mod-along-read plots.

## URLs consulted
- https://github.com/nanoporetech/modkit
- https://github.com/nanoporetech/modkit/blob/master/book/src/intro_summary.md
- https://api.github.com/repos/nanoporetech/modkit/releases/latest
- https://github.com/PacificBiosciences/pb-CpG-tools
- https://raw.githubusercontent.com/samtools/hts-specs/master/SAMtags.tex
- https://raw.githubusercontent.com/samtools/htslib/1.16/htslib/sam.h
- https://github.com/samtools/htslib/blob/develop/NEWS
