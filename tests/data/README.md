# SAMStat test fixtures

Hand-written inputs for known bugs and edge cases. Expected values below were
computed by hand from the SAM specification, NOT from current samstat output;
where current samstat differs, the fixture exposes a bug. Every SAM here parses
with htslib 1.16 `test/test_view` (the BAM files were made with it).

All SAM files have `@HD VN:1.6` and `@SQ SN:chr1 LN:2000` (`small.sam`: LN:5000).
Read sequences are the repeating pattern `ACGT...` unless stated. Qualities are
`I` (Phred 40) unless stated. Source line numbers refer to the tree at the time
of writing (branch maintenance/v2.3).

MAPQ bins used by the JSON expectations (stable keys, display labels are not
renamed in 2.3.0): `gt30` = MAPQ >= 30, `lt30` = 10..29, `lt10` = 0..9,
`unmapped` = flag 0x4. (Bin map: src/collect/collect.c:1214-1221.)

Note: macOS file systems are case-insensitive, so `READS.FQ` must not be
accompanied by a `reads.fq`; the lower-case FASTQ exists only inside `reads.fq.gz`.

## a) read_1024.sam

- Purpose: read of exactly 1024 bp (the initial conversion buffer size).
- Bug: src/htsinterface/htsglue.c:71,97 `b_len = 1024` and `if(s->len > b_len)`
  should be `>=`; `buffer[s->len] = 0` (line 109) then writes 1 byte past the
  1024-byte buffer.
- Expected: total reads 1; mapped 1; MAPQ 60 -> `gt30` = 1; read length 1024
  (one entry in the length distribution at 1024); mismatches 0 (MD:Z:1024);
  qualities all Q40 (1024 bases); no crash, exit 0.

## b) seq_star_softclip.sam

- Purpose: mapped record with `SEQ=*`, `QUAL=*`, CIGAR `10S40M`.
- Bug: src/collect/collect.c:~369 (and :~435, :~498 equivalents)
  `t[q_idx][itm->len - (clip_start + clip_end)]++` with len 0 and clip 10 gives
  index -10 (out-of-bounds write).
- Expected: total reads 1; MAPQ 40 -> `gt30` = 1; stored read length 0 (SEQ is
  absent, l_qseq = 0; a length of -10 or 40 is wrong); no mismatches; no
  quality values; no crash, exit 0.

## c) softclip5_md.sam

- Purpose: soft clip at start of CIGAR combined with MD mismatch.
- Bug: src/sambamparse/sam_bam_parse.c:135-136 `case 'S': sp++;` should advance
  by `Oplen`, so the read pointer is misaligned by 4 for every following base.
- Reads (both: SEQ `TTTTT ACGTACGTA C GTACGTACGT`, 25 bp, CIGAR `5S20M`,
  `NM:i:1 MD:Z:9A10`):
  - `sc5_fwd`, flag 0, POS 11.
  - `sc5_rev`, flag 16, POS 101.
- Expected: total reads 2 (1 per record), MAPQ 60 -> `gt30` = 2;
  read length 25; aligned length 20; exactly 1 mismatch, reference base `A`,
  read base `C`. Position: aligned offset 9 (0-based, i.e. 10th aligned base),
  SEQ index 14 (0-based), reference position 20 for `sc5_fwd` (11+9) and 110
  for `sc5_rev` (101+9). The mismatch is a single A->C substitution; with the bug
  the reported read base is taken from 4 positions to the left (`A`/`T`), so the
  mismatch appears as a match or as a wrong base.
- Clipping: `sc5_fwd` n_clip5 = 5, n_clip3 = 0. `sc5_rev` (reverse strand): after
  the strand swap at sam_bam_parse.c:254-256, n_clip5 = 0, n_clip3 = 5.
- Quality: 25 bases of Q40 each.

## d) hardclip_softclip.sam

- Purpose: hard clip adjacent to soft clip.
- Bug: src/sambamparse/sam_bam_parse.c:80-86 only sets n_clip5 when the `S` is
  CIGAR element 0 (and n_clip3 when it is the last element). With `5H10S40M`
  the `S` is element 1, so n_clip5 stays 0; likewise `40M10S5H`.
- Reads (SEQ 50 bp each, hard-clipped bases are not in SEQ):
  - `hs_5prime`, CIGAR `5H10S40M`: n_clip5 = 10, n_clip3 = 0.
  - `hs_3prime`, CIGAR `40M10S5H`: n_clip5 = 0, n_clip3 = 10.
- Expected: total reads 2; MAPQ 50 -> `gt30` = 2; read length 50 each; aligned
  length 40 each; 0 mismatches (no MD); quality Q40 x 100 bases.

## e) supplementary.sam

- Purpose: one read name with primary (flag 0, `50M`, MAPQ 60), supplementary
  (flag 2048, `30S20M`, MAPQ 60) and secondary (flag 256, `50M`, MAPQ 0).
- Bug: src/htsinterface/htsglue.c:76 only skips BAM_FSECONDARY; 0x800 records
  are counted as extra reads.
- Expected: total counted reads 1; `gt30` = 1; `lt10` = 0 (the MAPQ 0
  secondary is not counted); read length 50; supplementary alignments counted
  separately = 1; secondary alignments counted separately = 1 (once the
  counters exist); 0 mismatches; quality Q40 x 50.

## f) no_qual.sam

- Purpose: mapped reads with `QUAL=*`. Run with `--plotend --norm` (and also
  with default options).
- Bug: src/htsinterface/htsglue.c:117-120 fills quality with 0xFF; collect.c
  (:302 and similar) computes `q = 0xFF - 33 = 222` and clamps to 41, so these
  reads are reported as all-Q41.
- Reads: `nq_hi` MAPQ 40, `nq_lo` MAPQ 5, `nq_mid` MAPQ 20; each 30 bp, `30M`.
- Expected: total reads 3; `gt30` = 1, `lt30` = 1, `lt10` = 1; read length 30
  for all (3 reads at length 30); 0 mismatches; NO quality values recorded
  (quality histogram all zero; in particular Q41 count = 0); no NaN/inf in the
  report or JSON.

## g) high_qual.sam

- Purpose: qualities above Q41, which are clamped to the top bin.
- Reads (25 bp, `25M`, MAPQ 60, no mismatches):
  - `hq1`: QUAL = `I`x5 `J`x5 `K`x5 `]`x5 `~`x5, i.e. Q40 x5, Q41 x5, Q42 x5,
    Q60 x5, Q93 x5.
  - `hq2`: QUAL = `]` x25 (Q60).
- Bug site: src/collect/collect.c:302, 398, 469, 528 (`if(q > 41) q = 41`).
- Expected: total reads 2 (`gt30` = 2); 50 bases in all: 5 at Q40 and 45
  in the Q41 bin (`hq1`: 5 J + 5 K + 5 ] + 5 ~ = 20, `hq2`: 25). No out-of-range index, no
  crash. Per-position: `hq1` positions 0-4 Q40, 5-24 Q41(clamped); `hq2` 0-24
  Q41(clamped).

## h) md_inconsistent.sam and md_8digit.sam

- Purpose: MD tags that do not agree with the alignment, to test that MD
  parsing is bounded by the alignment length.
- Bug: src/sambamparse/sam_bam_parse.c:195-203 loops `c` times over
  `genome[pos]` with `pos` unbounded by aln_len (alloc: 256 minimum);
  :179 `char tmp_num[8]` overflows for a number of >= 8 digits.
- `md_inconsistent.sam` (all `20M`, MAPQ 60, Q40):
  - `md_too_long`: `MD:Z:1000` (claims 1000 matches for 20 aligned bases).
  - `md_mismatch_then_long`: `MD:Z:10A5000`; SEQ has `G` at aligned offset 10.
  - `md_ok_control`: `MD:Z:20`, read matches.
  - Expected: total reads 3 (`gt30` = 3); read length 20 each; no crash and no
    out-of-bounds access; the control read has 0 mismatches. Mismatches
    reported from the inconsistent reads must lie within the alignment
    (offset < 20); `md_mismatch_then_long` has exactly 1 mismatch at aligned
    offset 10 (ref `A`, read `G`); `md_too_long` has 0 mismatches. Nothing
    beyond offset 19 is reported.
- `md_8digit.sam`: one read, `20M`, `MD:Z:00000020` (valid per the SAM spec's
  `[0-9]+`; 8 digits). Expected: total reads 1; 0 mismatches; no crash.

## i) empty.sam and empty.fq

- `empty.sam`: header only, 0 records. `empty.fq`: zero-byte file.
- Expected: total reads 0; a report is written; a warning is printed; exit 0;
  no division by zero (no NaN/inf in report or JSON).

## j) reads.fq.gz, reads.fa.gz, READS.FQ

- Purpose: gzip-compressed (`gzip -n`) FASTQ/FASTA and an upper-case extension.
- Content (all three): 3 reads `read1`..`read3` of lengths 20, 25, 30
  (`ACGT` repeat); FASTQ quality all `I` (Q40).
- Expected: total reads 3; length distribution 20:1, 25:1, 30:1; all reads in
  the `unmapped` bin (sequence input has no mapping); FASTQ: Q40 x 75 bases;
  FASTA: no quality values. 75 bases in total. Exit 0.

## k) small.bam, truncated.bam

- `small.sam` -> `small.bam` (made with htslib `test_view -b`): 20 reads
  `small00`..`small19`, `40M`, positions 1,51,..., even index MAPQ 60, odd
  MAPQ 15, Q40, `MD:Z:40`. Expected: total reads 20; `gt30` = 10; `lt30` = 10;
  read length 40 (x20); 0 mismatches; Q40 x 800 bases; exit 0.
- `truncated.bam`: the first 200 of the 340 bytes of `small.bam`. The file is a
  single BGZF block, so the cut is inside that block (and loses the EOF marker);
  htslib reports `bgzf_read ... Read block operation failed`. Expected: samstat
  prints an error for this input, and exits 1 (processing other inputs still);
  it must not crash, hang or write a partial report that looks valid.

## l) not_a_bam.bam, unknown.xyz

- `not_a_bam.bam`: plain text with a `.bam` extension. Expected: clean error
  message ("unsupported or unknown format"-style), exit 1, no crash.
- `unknown.xyz`: a valid SAM (identical to `small.sam`) with an unknown
  extension. Expected: same as `small.sam` if content sniffing is used (20
  reads, exit 0); otherwise a clean, specific error about the unknown extension
  and exit 1. Either way, no crash. (The maintainer lead's decision on this is
  pending; tests should assert only no crash and exit code in {0,1}.)

## m) paired_small.sam

- 2 pairs: `pairA` flags 99/147 (MAPQ 60), `pairB` flags 99/147 (MAPQ 20);
  R1 reads are `30M`, R2 reads are `25M`; `MD:Z` all matches.
- Expected: total reads 4; read1 = 2; read2 = 2; `gt30` = 2 (pairA R1+R2),
  `lt30` = 2 (pairB R1+R2); R1 length 30 x2; R2 length 25 x2; 0 mismatches;
  Q40 x 110 bases (60 from R1, 50 from R2).

## n) softclip_shift.sam

- Purpose: soft clip whose length is not a multiple of 4, with a read that is not
  periodic. With `case 'S': sp++` the aligned bases are taken from the wrong place
  in SEQ; `softclip5_md.sam` clips 5 and its sequence repeats every 4 bases, so
  it cannot show that.
- Read `sc3_shift` (`3S20M`, SEQ `GGG` + `ACGT` x5, `NM:i:1 MD:Z:2T17`).
- Expected: total reads 1; length 20; exactly 1 mismatch at aligned offset 2:
  reference `T`, read base `G` (SEQ index 5). With the bug the read base is `A`.

## o) md_long_number.sam

- Purpose: MD digit runs far longer than any integer (50 and 34 digits).
- `md_long_digits` (`20M`, `MD:Z:99...9`): claims an absurd number of matches.
  `md_mm_long_digits` (`20M`, `MD:Z:10A99...9`): one mismatch at aligned offset
  10 (ref `A`, read `G`), then an absurd number of matches.
- Expected: total reads 2 (`gt30` = 2); read length 20; the first read has 0
  mismatches, the second exactly 1 (read base `G` at offset 10); no crash and no
  out-of-bounds access (visible with `-DSAMSTAT_SANITIZE=ON`).

## p) zero_bytes.sam

- Purpose: a zero byte file with a `.sam` extension (no header at all, unlike
  `empty.sam`).
- Expected: total reads 0; a report is written; exit 0. It is the only input
  without a readable header that is accepted (`not_a_bam.bam` is rejected).
