#ifndef PRECHECK_H
#define PRECHECK_H

/* Streaming pre-check of text inputs (FASTA, FASTQ, SAM; plain, gzip or bgzip).

   htslib reads a FASTA/FASTQ record, and a SAM or FASTQ line, into memory in one piece
   before samstat sees it. The read length limit (reader_opts.max_read_len) is applied
   after sam_read1() returns, so a genome FASTA costs about a chromosome of memory per
   record before it is rejected, and a compressed file with one giant line and no newline
   grows memory without bound. precheck_input() reads the start of the input a second time,
   in fixed-size chunks (constant memory), and rejects it before htslib loads a record:

   - FASTA/FASTQ: a record with more than max_read_len bases (the same rule as the reader,
     so nothing the reader accepts is rejected here), or a line longer than
     max_read_len + PRECHECK_LINE_SLACK (a header line; a sequence or quality line that long
     belongs to a record that is too long or is malformed).
   - SAM: a header line longer than max_read_len + PRECHECK_LINE_SLACK, or a record whose
     SEQ field has more than max_read_len bases. A SAM record line itself is not limited:
     SEQ, QUAL and tags of a read at the limit legitimately take several times its length.

   Only the first PRECHECK_SCAN_MAX bytes of decompressed content are checked: a record or
   line that starts later is still found by the reader's own check, after htslib has
   loaded it. Standard input ("-") and other inputs that are not regular files (a pipe,
   process substitution) cannot be read twice and are not checked. BAM and CRAM are not
   checked (their records carry their length; htslib does not read them line by line). */

#include "hts.h"
#include <stdint.h>

#ifdef PRECHECK_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

/* bytes of decompressed content scanned at most */
#define PRECHECK_SCAN_MAX (64 * 1024 * 1024)
/* a line may be this much longer than max_read_len (room for a header line) */
#define PRECHECK_LINE_SLACK 65536

/* <format> is the format htslib detected for <filename>. Returns OK if the input passes or
   is not checked (wrong format, not a regular file, cannot be read: the reader reports
   that), FAIL with an error message if it is rejected. */
EXTERN int precheck_input(const char *filename, enum htsExactFormat format, uint64_t max_read_len);

#undef PRECHECK_IMPORT
#undef EXTERN

#endif
