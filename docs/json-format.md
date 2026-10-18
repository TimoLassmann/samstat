# SAMStat JSON report (format_version 2)

`samstat --json file.bam` writes `file.bam.samstat.json` next to the HTML report
(in the directory given with `-d` if used). The file holds the integer counts SAMStat
collected on fixed axes, not the data of the plots as drawn. Formal description:
[`samstat.schema.json`](samstat.schema.json) (JSON Schema draft 2020-12). The layout is
frozen for all 2.x releases.

Properties of the file:

* Deterministic: same input and version give byte-identical output apart from
  `meta`. No timestamps, no absolute paths. Key order is fixed.
* Every key is always present (no existence checks needed), except `read2` (see below)
  and, in files written by 2.4.0, `mismatch_base` (added later).
* All arrays hold exact non-negative integers and are trimmed after their last non-zero
  entry: an index past the end of an array means 0, and a class without data gives `[]`.
* Floating point values appear only in `summary`, are printed with four decimals and are
  computed from integer accumulators. A value that does not apply is `null`.
* The size does not depend on the read length: every array is bounded (below). A
  1 Mb read gives a file well under 1 MB.
* Format changes that break consumers bump `format_version`; new optional fields do not.

### Changes from format_version 1 (2.3)

Consumers of version 1 must change:

* `format_version` is `2` and the layout is new. The transitional `v2` key of the 2.4
  development builds is gone; its members are the top level.
* Removed: `read*.base_composition`, `read*.base_quality` (the Q x position matrix),
  `read*.mismatches` (per-base mismatches), `read*.insertions`, `read*.deletions`, and the
  clip-trimmed `read*.length`. Their replacements are in `read*.absolute` / `read*.relative`
  (composition, mismatch, insertion and deletion counts per position bin) and `read*.length`
  (now the full read length, clipped bases included).
* New `summary` members, new `axes`, and floating point values in `summary`.
* `summary.is_aligned` is `true` if the header has an `@SQ` line or at least one primary
  read is mapped (v1: true for every SAM/BAM). It is always `false` for FASTA/FASTQ.

## Top level

| key | type | description |
|---|---|---|
| `format_version` | integer | Always `2` for this layout. |
| `meta` | object | The only part that may differ between runs. |
| `summary` | object | File-level numbers. |
| `mapq_bins` | array | The four mapping-quality classes, in fixed order, with stable ids. |
| `axes` | object | What the indexes of the arrays mean. |
| `read1` | object | Statistics for read 1 (all reads of unpaired data). |
| `read2` | object | Same layout as `read1`. Present only if the input has read 2 records (`summary.read2_reads > 0`). |

### `meta`

| key | description |
|---|---|
| `samstat_version` | SAMStat version string, e.g. `"2.4.0"`. |
| `input_file` | Base name of the input file (no directory); `"stdin"` for standard input. |

### `summary`

All statistics are over primary records; secondary and supplementary alignments are
counted separately and are not reads.

| key | description |
|---|---|
| `total_reads` | Reads counted (read 1 + read 2). |
| `read1_reads`, `read2_reads` | Reads counted as read 1 / read 2. Records without a read 2 flag, and all FASTA/FASTQ records, count as read 1. |
| `paired` | Reads flagged as paired in sequencing (0x1). |
| `properly_paired` | Reads flagged as properly paired (0x2). |
| `secondary`, `supplementary` | Secondary (0x100) and supplementary (0x800) alignments seen. With `-p` only those up to the cutoff. |
| `is_aligned` | `true` if the header has at least one `@SQ` line or at least one primary record is mapped; always `false` for FASTA/FASTQ. When `false` every read is in `unmapped` and the alignment statistics are empty. |
| `is_partial_report` | `true` if reading stopped early because of `-p/--peek` and at least one more primary read follows the cutoff (or the rest of the file cannot be read); counts then cover only the reads seen. To find out, samstat reads on past the cutoff (through any secondary and supplementary records, without counting them, up to the end of the file if only those follow); `false` when no primary read follows. |
| `total_bases`, `min_length`, `max_length` | Bases and shortest / longest read, full read length (clipped bases included). `min_length` is 0 when there are no reads. |
| `n50`, `median_length` | Lower edge of the storage length bin where the cumulative fraction of bases (N50, counted from the longest reads) or of reads (median) first reaches 50%. Exact below 1000 bp. |
| `mean_q` | Mean base quality over all bases with a quality (clipped bases included), taken over error probabilities, not over Phred values. `null` without qualities. |
| `pct_mapped`, `pct_mapq_ge30` | Percent of reads that are mapped / mapped with MAPQ >= 30. `null` for input that is not aligned. |
| `pct_clipped` | Percent of the bases of mapped reads that are soft-clipped. `null` without mapped reads. |
| `median_identity` | Lower edge of the median 0.1% bin of the gap-compressed identity, in percent. `null` without a computable identity. |
| `layout` | `"long"` if `median_length >= 1000`, else `"short"`. Decides which position axis the HTML report draws; the JSON carries both. |
| `reads_without_md` | Mapped reads with neither an `MD` tag nor `X` operations in the CIGAR (they contribute nothing to `aligned_md` and `mismatches`). |
| `identity_reads`, `identity_na` | Reads in the identity histogram; mapped reads whose identity cannot be computed (no `NM` tag and no `=`/`X` operations, or an inconsistent `NM`). |
| `bases_beyond_position_axis` | Aligned bases at read positions >= 1000, which the absolute axis does not hold. |
| `mapq_counts` | `{"read1": {...}, "read2": {...}}`, reads per mapq id. `read2` only if present. The counts of both mates sum to `total_reads`. |

Identity is `(M - X) / (M + I + D)` where `M` counts the bases of `M`, `=` and `X`
operations, `I` and `D` the insertion and deletion events (the `de` convention of minimap2)
and `X` the mismatches from `NM` (or from `X` operations if `NM` is absent).

## Mapq bins

`mapq_bins` is an array of `{"id": ..., "label": ...}` in this order. Use the `id`
as a key; the `label` is display text and may change in a later release.

| id | contents |
|---|---|
| `mapq_ge30` | mapped, MAPQ >= 30 (MAPQ 255 included) |
| `mapq_10_29` | mapped, MAPQ 10 to 29 |
| `mapq_lt10` | mapped, MAPQ < 10 |
| `unmapped` | unmapped records, and all FASTA/FASTQ records |

All per-class objects below have exactly these four ids as keys.

## `axes`

| key | description |
|---|---|
| `length_edges` | Lower edge of each of the 1922 storage length bins: lengths 0 to 999 have a bin each, from 1000 bp the edges grow by about 1% per bin (`e -> max(e + 1, (101 e + 99) / 100)`) up to the first edge >= 10,000,000. That last bin holds every longer read. |
| `absolute` | `{"bins": 1000, "label": "position in read (bp)"}`: bin `i` is read position `i` (0-based, read orientation, clipped bases excluded). Positions >= 1000 are not counted. |
| `relative` | `{"bins": 100, "label": "position in read (% of length)"}`: bin `(100 * i) / L` for position `i` of an aligned part of length `L`. |
| `identity_bins` | 1001: bin `b` holds identity in [b/10 %, (b+1)/10 %), bin 1000 exactly 100%. |
| `quality_error_scale` | 2^30 = 1073741824, the scale of `qual_error_sum`. |

## `read1` / `read2`

| key | shape | description |
|---|---|---|
| `n_reads` | integer | Reads of this read number. |
| `length` | `{"reads": {id: [int]}, "bases": {id: [int]}}` | Reads, and the sum of their lengths, per storage length bin (index = bin of `axes.length_edges`). The length is the full read length. |
| `read_mean_q` | `{id: [int]}` | Reads per mean quality; index = Q 0 to 93, bin `q` = read mean Q in [q, q+1). Reads without qualities are not counted. At most 94 entries. |
| `identity` | `{id: [int]}` | Mapped reads per identity bin (0.1% steps). `[]` for `unmapped`. At most 1001 entries. |
| `absolute` | position block | Counts per read position, up to 1000 bins. |
| `relative` | position block | Counts per percent of the aligned read length, up to 100 bins. |

A **position block** holds, per mapq id and position bin:

| key | description |
|---|---|
| `composition` | `{id: {"A": [int], "C": [int], "G": [int], "T": [int], "N": [int]}}`: read bases; anything that is not A, C, G or T counts as N. |
| `qual_bases` | Bases with a quality value. |
| `qual_error_sum` | Sum of the error probabilities of those bases, in units of 1 / `axes.quality_error_scale` (`E[q] = round(10^(-q/10) * 2^30)`, Phred values above 93 count as 93). Mean Q of a bin is `-10 log10(qual_error_sum / qual_bases / 2^30)`. |
| `aligned` | Aligned read bases (CIGAR `M`, `=`, `X`, `I`) of mapped reads: denominator for insertions and deletions. |
| `aligned_md` | The same, only for reads with an `MD` tag or `=`/`X` operations: denominator for mismatches. |
| `mismatches` | Mismatching bases (from `MD`, or from `X` operations when `MD` is absent; a read of `=` operations only has none). |
| `mismatch_base` | `{id: {"A": [int], "C": [int], "G": [int], "T": [int], "N": [int]}}`: the same mismatches split by the **read** base at that position (read orientation, so a reverse-strand read counts the complement of its stored base; anything that is not A, C, G or T counts as N). The five arrays add up to `mismatches` in every bin; `unmapped` has empty arrays. It works without `MD` (CIGAR `X`). Added after 2.4.0 as a new optional field (`format_version` stays 2): files written by 2.4.0 do not have it. |
| `insertions` | Insertion events starting at this read position (one per gap opening). |
| `deletions` | Deletion events at this read position (one per gap opening). |

On aligned input the position axes cover the aligned part of the read (soft-clipped bases
are excluded); on FASTA/FASTQ they cover the whole read.

Array length bounds: absolute-axis arrays at most 1000, relative-axis arrays at most 100,
identity arrays at most 1001, `read_mean_q` at most 94.

## Example (abridged, `dev/xx#pair.sam`)

```json
{
"format_version": 2,
"meta": {"samstat_version": "2.4.0", "input_file": "xx#pair.sam"},
"summary": {
"total_reads": 6,
"...": "...",
"mean_q": 9.0000,
"median_identity": null,
"layout": "short",
"mapq_counts": {"read1": {"mapq_ge30": 0, "mapq_10_29": 0, "mapq_lt10": 3, "unmapped": 0}, "read2": {"...": 0}}
},
"mapq_bins": [ {"id": "mapq_ge30", "label": "MAPQ >= 30"}, "..." ],
"axes": {"length_edges": [0,1,2,"..."], "absolute": {"bins": 1000, "label": "position in read (bp)"}, "...": "..."},
"read1": {
"n_reads": 3,
"length": {"reads": {"mapq_ge30": [], "mapq_lt10": [0,0,0,0,0,0,0,0,0,0,3], "...": []}, "bases": {"...": []}},
"read_mean_q": {"mapq_lt10": [0,0,0,0,0,0,0,0,0,3], "...": []},
"...": "..."
}
}
```

(The real file is not abridged and has one object per line.)

## Testing and regeneration

`tests/cmake/check_json.cmake` checks individual values by path
(`read1.length.reads.mapq_lt10.10=3`, `mapq_bins#len=4`), and `tests/cmake/compare_golden.cmake`
compares a whole report with a golden file in `tests/golden/expected/` (ignoring `meta`).
After an intended format change regenerate the goldens with the `update-goldens` target;
see the header of `tests/json/CMakeLists.txt`. `tests/cmake/check_schema_keys.cmake`
compares the member names of a report with `samstat.schema.json`.
