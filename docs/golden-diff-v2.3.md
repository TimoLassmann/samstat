# Golden report changes in 2.3.0

The golden reports in `tests/golden/expected/` were captured before the 2.3.0
fixes (commit `198270d`, "baseline goldens (pre-fix behaviour)") and regenerated
after all work packages were merged (`cmake --build <build> --target update-goldens`).
This page explains every difference. Each change is traced to an intended fix.

## Changed in every golden file

| Field | Change | Fix |
| --- | --- | --- |
| `summary.secondary` | added (new key) | Secondary (0x100) and supplementary (0x800) records are counted separately (wp/collect, wp/reader) |
| `summary.supplementary` | added (new key) | same |

For all inputs other than those listed below, these two keys are the only
change (`summary.secondary` and `summary.supplementary` are 0). The following
files showed no other difference: `aux#aux`, `aux.fa`, `c1#*`, `c1.fa`,
`ce#1`, `ce#2`, `ce#5b`, `ce#tag_*`, `ce#unmap*`, `empty`, `fieldarith`,
`high_qual`, `md_inconsistent`, `no_qual`, `paired_small`, `read_1024`,
`small.bam`, `small.sam`, `softclip5_md`, `xx#*`, `xx.fa`.

`no_qual`, `high_qual`, `md_inconsistent`, `read_1024` and `softclip5_md` did
not change in the JSON. Their fixes are still covered by dedicated tests
(missing quality `*` and clamping affect only the `--plotend`/`--norm` HTML
plots, the MD overflow and the 1024 base buffer overflow only show up as
crashes or sanitizer reports, and `softclip5_md` has a 5' clip that is a
multiple of the sequence period).

## Files with further changes

### `ce#5.sam.json` — supplementary exclusion

- `summary.total_reads` 6 -> 5, `summary.read1_reads` 6 -> 5,
  `summary.mapq_counts.read1.mapq_lt10` 6 -> 5, `read1.n_reads` 6 -> 5.
- `summary.supplementary` 0 -> 1.
- `read1.length.mapq_lt10` bin 100: 6 -> 5, and base composition, base quality
  and deletion arrays of `mapq_lt10` lose the contribution of that one record.
- Fix: supplementary alignments (0x800) are no longer counted as reads
  (wp/reader). The input contains one such record.

### `supplementary.sam.json` — supplementary exclusion and soft-clip length

- `summary.total_reads` 2 -> 1, `read1_reads` 2 -> 1,
  `mapq_counts.read1.mapq_ge30` 2 -> 1, `summary.secondary` = 1,
  `summary.supplementary` = 1.
- `read1.length.mapq_ge30`: the 20-base entry (the `30S20M` supplementary
  record) is gone; one 50-base read remains.
- Base composition and quality arrays lose that record's bases.
- Fix: supplementary exclusion. The secondary record (0x100) was already
  excluded, it is now counted in `summary.secondary`.

### `hardclip_softclip.sam.json` — soft clip after hard clip

- `read1.length.mapq_ge30`: both reads move from length 50 to length 40.
- `read1.base_composition.mapq_ge30` and `read1.base_quality.mapq_ge30` shift
  by the clip length: per-position composition now starts at the first aligned
  base, not at the first base of the stored sequence.
- Fix: for `5H10S40M` and `40M10S5H` the soft clip is located past the `H`,
  so the 10 soft-clipped bases are removed from the aligned length and from
  the per-position statistics (wp/reader, "set clips around H").

### `seq_star_softclip.sam.json` — SEQ `*` with a soft-clip CIGAR

- `read1.length.mapq_ge30`: before, no length was recorded (`[]`, the index was
  negative, an out-of-bounds write); now bin 0 holds 1 read.
- Fix: records with SEQ `*` have an empty sequence and aligned length 0;
  length-histogram indexing is clamped (wp/reader and wp/collect).

## New golden files (no previous baseline)

`md_8digit.sam.json`, `md_long_number.sam.json`, `softclip_shift.sam.json` and
`zero_bytes.sam.json` belong to fixtures added by this release. They have no
pre-fix baseline; they record the post-fix output.

## Removed golden file

`not_a_bam.bam.json` is removed. The input is not a valid BAM; since 2.3.0
samstat reports an error, writes no report and exits 1. This is covered by
`reader:error:not_a_bam`.

## Unexplained changes

None. Every changed field above maps to one of the intended fixes.
