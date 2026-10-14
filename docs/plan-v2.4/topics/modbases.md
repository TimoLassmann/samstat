# Topic: modified bases (MM/ML) in SAMStat 2.4

Status: FINAL design. It is based on the "minimal" proposal, with a few ideas taken from the "researcher" proposal: the flag semantics, the captions, the empty-MM ASAN test, the uBAM test and schema validation both ways. It also corrects one claim in the minimal proposal (see section 4).

## 1. Decision

**Rich modbase QC stays OUT of 2.4. One conditional two-row check is IN.**

- There is no new plot, no CLI option and no new section.
- If at least one primary record carries an `MM` tag, the existing summary table gets two extra rows:
  1. **Reads with base-modification calls (MM/ML):** the count and the percentage of primary reads.
  2. **Reads with inconsistent MM/ML:** a count. It shows CHECK when it is above zero.
- If no record carries `MM`, the HTML, JSON and `_mqc` output are byte-identical to 2.3. The goldens do not change.

### Why this cut

- The modBAM failure a general QC tool should catch is silent loss: tags that were dropped, or tags that went stale after SEQ was trimmed or hard-clipped. Two counts detect it at almost no cost. Source: `research/modbases.md`, section 5.
- Call counts, fraction modified and probability histograms already exist in modkit (`modbam summary`, `sample-probs`). See https://github.com/nanoporetech/modkit and https://github.com/nanoporetech/modkit/blob/master/book/src/intro_summary.md. Copying them would add rows that are biology, not QC, and that maintenance falls on one maintainer.
- SAMStat's differentiator is error profiles by MAPQ, and modbase rows do not strengthen it. Nobody has asked for this: none of SAMStat's open issues mention modifications (`research/user-needs.md`).

## 2. Implementation: no htslib modification API, no htslib bump

The check uses only `bam_aux_get` and plain string scans, so it works on the bundled htslib 1.16 and on any later version.

- **Verified in the code:** `bam_parse_basemod` and `bam_next_basemod` exist in the bundled 1.16 headers (`build-asan/htslib/src/htslib/htslib/sam.h:2207`). This design deliberately does not use them.
- **Why not:** they need MN support (htslib 1.18), and an empty-MM segfault was fixed only in 1.23 (https://github.com/samtools/htslib/blob/develop/NEWS).
- **Unverified:** whether 1.16 crashes on an empty MM. The question is moot here, because this code never calls the htslib parser.
- **Bump policy:** this feature must not depend on the htslib bump that the perf topic proposes, and it must not justify one.

### Where it goes

Verified in the code: `read_bam_chunk` in `src/htsinterface/htsglue.c:82-89` already separates secondary and supplementary records (it counts them into `f_handle->n_secondary` and `f_handle->n_supplementary`) from primary records, and it handles primaries while the `bam1_t` is still available. `samstat.c:220` copies those counters into the stats struct, and `stat_report.c:473` renders them as conditional rows.

The new counters follow the same path:

- `f_handle->n_mm_reads` and `f_handle->n_mm_bad` are filled in the primary branch of `read_bam_chunk`.
- They are copied into `struct stat_collection` next to `n_secondary`. With the per-thread accumulators planned by the perf topic, they merge by addition.

### Per primary record

1. `const uint8_t *mm = bam_aux_get(b, "MM");`. If that returns NULL, do nothing. Records that only have the legacy `Mm`/`Ml` spellings are ignored (see DECISION NEEDED 2).
2. If `MM` is present, increment `n_mm_reads`.
3. Mark the read inconsistent and increment `n_mm_bad` if any of the following holds:
   - a. **MN mismatch.** `MN` is present and is not equal to `b->core.l_qseq`. MN is defined in the SAM tags spec, https://raw.githubusercontent.com/samtools/hts-specs/master/SAMtags.tex. The date it was introduced (2024) is UNVERIFIED.
   - b. **ML count mismatch.** `ML` is missing although MM lists at least one call, or the length of the ML array differs from the total number of comma-separated deltas in MM. A segment with no deltas, such as `C+m;` or `C+m?;`, contributes zero, so an empty MM with an absent or empty ML is consistent.
   - c. **Positions beyond SEQ.** In each segment, the sum of (delta + 1) is greater than the number of that canonical base in SEQ, read in its original orientation. When `BAM_FREVERSE` is set, count the complement in the stored SEQ. A canonical base of `N` matches every base.
4. Records with a `*` SEQ (`l_qseq == 0`) and an MM tag count as inconsistent.

**Correction to the minimal proposal.** Check (c) is the essential one. Hard-clipping or trimming SEQ leaves the MM and ML strings unchanged, so check (b) cannot see it. Without MN, only check (c) catches stale tags. Check (b) catches broken writers. The minimal proposal said (c) could be cut first. It must stay.

**Cost.** Reads without MM pay one aux lookup. Reads with MM pay one pass over the MM string and one pass over SEQ. That pass over SEQ is O(read length), but SAMStat already makes 10 to 14 passes per base (`research/code-collect.md`). Not benchmarked (UNVERIFIED), but expected to be negligible.

**Scope.** Primary records only, which includes unaligned uBAM records. Secondary and supplementary records are skipped, so there is no double counting. FASTA and FASTQ input never produce these rows.

## 3. Output

### HTML

The two rows sit in the existing summary or flags table and appear only when `n_mm_reads > 0`.

| Row | Value | Flag | Caption (one line) |
|---|---|---|---|
| Reads with base-modification calls | count and % of primary reads | info, never flagged | "Reads carrying MM/ML modification tags. 0% is normal if the basecaller was not run with a modification model." |
| Reads with inconsistent MM/ML | count | CHECK if > 0, otherwise OK | "Should be 0. Non-zero usually means the sequence was trimmed or hard-clipped after basecalling. Filter or regenerate the tags before running modkit." |

- CHECK uses the OK / CHECK / FAIL convention: a shape and a word, never colour alone. It is never FAIL, because the flag is advisory (`research/presentation.md`).

### JSON

- Add one optional object, `"modifications": {"reads_with_mm": N, "reads_inconsistent": M}`. It is present only when `reads_with_mm > 0`.
- The change is additive: list it as an optional property in `docs/samstat.schema.json` and document it in `docs/json-format.md`. Its size is constant.

### MultiQC `_mqc`

- Add two general-stats columns, `mm_read_pct` and `mm_inconsistent_reads`, written only for samples that have them.
- Samples without MM produce the same files as in 2.3. Document this in `docs/multiqc.md`.

## 4. Explicitly OUT for 2.4

- Per (base, mod code) call counts, % modified at p ≥ 0.5, the ML probability histograms, and modification by MAPQ bin. These would only be reconsidered under DECISION NEEDED 1.
- Adaptive or user-set thresholds, per-site pileups and bedMethyl output, motif or region filters, haplotype splits, and modification-along-read plots. Use modkit (https://github.com/nanoporetech/modkit) or pb-CpG-tools (https://github.com/PacificBiosciences/pb-CpG-tools) for these.
- Any `--mods` flag. The presence of the tag is the trigger.
- Use of the htslib modification API, and any htslib bump made for this feature.

## 5. Affected files

| File | Change |
|---|---|
| `src/htsinterface/htsglue.c` (+ header with `struct sam_bam_file`) | ~60 lines: the tag checks in the primary branch, two counters, and resetting them where `n_secondary` is reset (`htsglue.c:224`) |
| `src/collect/collect.h`, `collect.c` | two `uint64_t` fields with init and merge |
| `src/samstat.c` | copy the counters into the stats struct (next to line 220) |
| `src/report/stat_report.c` | the two conditional rows (pattern at line 473) |
| `src/report/json_report.c`, `docs/json-format.md`, `docs/samstat.schema.json` | the optional `modifications` object |
| `src/report/mqc_report.c`, `docs/multiqc.md` | the conditional columns |

Effort is about 1 day: roughly 120 lines of C, with about half the time going to fixtures and tests.

## 6. Tests

New fixtures go in `tests/data/longread/`. They are SAM text with hand-built MM/ML/MN tags, written by the seeded stdlib generator (`uv run python`) from the code-tests topic and converted to BAM with the bundled htslib if needed. Expected values come from the generator's ground truth, never from SAMStat output.

1. **Valid modBAM** (`C+m` and `C+h` on CpG positions, some reverse-strand reads): `reads_with_mm` equals the number of primaries, and `reads_inconsistent` is 0. This guards check (c) against a reverse-strand false positive.
2. **Mixed file:** half the reads carry MM, so `mm_read_pct` is 50.
3. **Stale tags.** Each case gets its own read, and the total must equal the generator's count.
   - MN ≠ l_qseq.
   - ML count ≠ MM count.
   - A primary whose SEQ was truncated with MM unchanged and no MN. Only check (c) can catch this one.
4. **Empty MM** (`C+m;`, `C+m?;`) with an empty or absent ML: counted as present and consistent. Run it in the ASAN job.
5. **Secondary and supplementary records carrying MM** are excluded from both counts.
6. **Unaligned uBAM with MM** is counted.
7. **No-MM regression:** all existing goldens are unchanged. `check_json` asserts that the `modifications` key is absent on non-MM input and present on modBAM input. The JSON validates against the schema in both cases.
8. **Report presence:** the row labels appear in the HTML only for the modBAM fixtures.
9. **`-t 1` vs `-t 4`:** identical output, which tests the additive merge once the per-thread accumulators exist.
10. **MultiQC:** the columns appear only for the modBAM sample.

## 7. Revisit trigger

Reconsider per-code call counts and % modified at p ≥ 0.5 only when both of these hold:

- htslib is at 1.23 or later for other reasons, so `bam_next_basemod` is safe to use.
- A user asks for modification QC in SAMStat.

Adding them later is additive: a new JSON key and new rows.

## 8. DECISION NEEDED

1. **Per-code % modified.** The researcher proposal's item 3 adds per (base, code) calls and % modified at p ≥ 0.5, with a cap of 16 codes. Is it in 2.4 if the perf topic's htslib 1.24 bump lands anyway?
   - Recommendation: no. Keep 2.4 at two rows.
2. **Legacy `Mm`/`Ml` tags.** Count them as MM, or ignore them? How common they still are is UNVERIFIED.
   - Recommendation: ignore them.
3. **Flag level.** Is CHECK on any non-zero inconsistent count right, or should there be a small tolerance, for example 0.1% of MM reads?
   - Recommendation: flag any non-zero count. Valid writers produce 0.

## 9. Unverified

- Whether htslib 1.16 crashes on an empty MM. This design does not depend on the answer.
- When MN was introduced to the spec.
- How common `Mm`/`Ml` tags still are.
- Runtime cost of the scan, which has not been benchmarked.
- Whether real Dorado or Revio pipelines produce stale MM on primary records often enough for this check to fire. Check this against a trimmed real modBAM before release.
