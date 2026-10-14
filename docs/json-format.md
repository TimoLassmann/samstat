# SAMStat JSON report (format_version 1)

`samstat --json file.bam` writes `file.bam.samstat.json` next to the HTML report
(in the directory given with `-d` if used). The file holds the raw counts SAMStat
collected, not the interval-binned data shown in the plots. Formal description:
[`samstat.schema.json`](samstat.schema.json) (JSON Schema draft 2020-12).

Properties of the file:

* Deterministic: same input and version give byte-identical output apart from
  `meta`. No timestamps, no absolute paths. Key order is fixed.
* All counts are exact non-negative integers. There are no floating point values.
* Format changes that break consumers bump `format_version`; new optional fields do not.

## Top level

| key | type | description |
|---|---|---|
| `format_version` | integer | Always `1` for this layout. |
| `meta` | object | The only part that may differ between runs. |
| `summary` | object | File-level counts. |
| `mapq_bins` | array | The four mapping-quality bins, in fixed order, with stable ids. |
| `read1` | object | Per-position statistics for read 1 (all reads of unpaired data). |
| `read2` | object | Same layout as `read1`. Present only if the input has read 2 records (`summary.read2_reads > 0`). |

### `meta`

| key | description |
|---|---|
| `samstat_version` | SAMStat version string, e.g. `"2.3.0"`. |
| `input_file` | Base name of the input file (no directory). |

### `summary`

| key | description |
|---|---|
| `total_reads` | Reads counted (read 1 + read 2). |
| `read1_reads`, `read2_reads` | Reads counted as read 1 / read 2. Records without a read 2 flag, and all FASTA/FASTQ records, count as read 1. |
| `paired` | Reads flagged as paired in sequencing (0x1). |
| `properly_paired` | Reads flagged as properly paired (0x2). |
| `secondary` | Secondary alignments (0x100) seen. They are excluded from all per-read statistics and `total_reads`. |
| `supplementary` | Supplementary alignments (0x800) seen. Excluded from per-read statistics and `total_reads`, like secondary alignments. |
| `is_aligned` | `true` for SAM/BAM input, `false` for FASTA/FASTQ. It reflects the file type only, so it is also `true` for a SAM/BAM without reads; check `total_reads` for that. When `false` the mapq bins carry no meaning: every read is in `unmapped`, and `mismatches`, `insertions` and `deletions` are empty. |
| `is_partial_report` | `true` if reading stopped early because of `-p/--peek`; counts then cover only the reads seen. |
| `mapq_counts` | `{"read1": {...}, "read2": {...}}`, reads per mapq bin (keys below). `read2` only if present. |

Within one read, a bin's count is the number of reads that fall in that bin; the bin
counts of `read1` and `read2` sum to `total_reads`.

## Mapq bins

`mapq_bins` is an array of `{"id": ..., "label": ...}` in this order. Use the `id`
as a key; the `label` is the display text used in the HTML report and may change in
a later release.

| id | label (2.3.0) | contents |
|---|---|---|
| `mapq_ge30` | `MAPQ gt 30` | mapped, MAPQ >= 30 |
| `mapq_10_29` | `MAPQ lt 30` | mapped, MAPQ 10 to 29 |
| `mapq_lt10` | `MAPQ lt 10` | mapped, MAPQ < 10 |
| `unmapped` | `Unmapped` | unmapped records, and all FASTA/FASTQ records |

All per-bin objects below have exactly these four ids as keys.

## `read1` / `read2`

| key | shape | description |
|---|---|---|
| `n_reads` | integer | Reads of this read number. |
| `length` | `{bin: [int]}` | Length histogram: `length[bin][L]` is the number of reads with length `L` (aligned length excluding soft/hard clipped bases for SAM/BAM). |
| `base_composition` | `{bin: {A,C,G,T,N: [int]}}` | `base_composition[bin][base][i]`: count of `base` at position `i` (0-based, clipped bases excluded). Anything that is not A, C, G or T is counted as `N`. |
| `base_quality` | `{bin: [[int]]}` | `base_quality[bin][q][i]`: count of bases with Phred quality `q` (0 to 41, higher values are clipped to 41) at position `i`. See below for the shape. |
| `mismatches` | `{bin: {A,C,G,T,N: [int]}}` | `mismatches[bin][base][i]`: reads with a mismatch at read position `i`, keyed by the read base. Needs an `MD` tag in SAM/BAM. |
| `insertions` | `{bin: [int]}` | `insertions[bin][i]`: insertions starting at read position `i`. |
| `deletions` | `{bin: [int]}` | `deletions[bin][i]`: deletions at read position `i` (positions are in read coordinates). |

### Trimmed arrays

Arrays are trimmed after their last non-zero entry; entries past the end are zero.
Rows of the same bin and section share one length, so for `base_composition`,
`base_quality` and `mismatches` all bases (or all Q values) of a bin have the same
number of positions. A bin with no data yields empty arrays (`[]`) instead of arrays
of zeros. Consequences for consumers:

* An index past the end of an array means 0.
* `base_quality[bin]` is either `[]` (no quality data: FASTA input, or an empty bin)
  or exactly 42 arrays, one per Q value 0 to 41.
* The keys themselves are always present, so no existence checks are needed.

Trimming keeps the file small for the common case, but a single very long read
dominates the size: the quality matrix has 42 rows per position, so one read of
1 Mb produces a file of around 100 MB.

## Example (abridged, `dev/xx#pair.sam`)

```json
{
"format_version": 1,
"meta": {"samstat_version": "2.3.0", "input_file": "xx#pair.sam"},
"summary": {
"total_reads": 6, "read1_reads": 3, "read2_reads": 3,
"paired": 6, "properly_paired": 6,
"is_aligned": true, "is_partial_report": false,
"mapq_counts": {
"read1": {"mapq_ge30": 0, "mapq_10_29": 0, "mapq_lt10": 3, "unmapped": 0},
"read2": {"mapq_ge30": 0, "mapq_10_29": 0, "mapq_lt10": 3, "unmapped": 0}
}
},
"mapq_bins": [ {"id": "mapq_ge30", "label": "MAPQ gt 30"}, "..." ],
"read1": {
"n_reads": 3,
"length": {"mapq_ge30": [], "mapq_10_29": [], "mapq_lt10": [0,0,0,0,0,0,0,0,0,0,3], "unmapped": []},
"...": "..."
}
}
```

(The real file is not abridged and has one object per line, as above.)

## Testing and regeneration

`tests/cmake/check_json.cmake` checks individual values by path
(`read1.length.mapq_lt10.10=3`, `mapq_bins#len=4`), and `tests/cmake/compare_golden.cmake`
compares a whole report with a golden file in `tests/json/golden/` (ignoring `meta`).
After an intended format change regenerate a golden with `-DUPDATE=ON`; see the header
of `tests/json/CMakeLists.txt`.
