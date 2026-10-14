#ifndef HTSGLUE_H
#define HTSGLUE_H

#include "tld.h"
#include "sam.h"
#include "hts.h"
#include <htslib/thread_pool.h>
#include <stdint.h>

#ifdef HTSGLUE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct aln_data {
        tld_strbuf* md;
        uint32_t* cigar;
        uint8_t* genome;
        uint8_t* read;
        /* uint8_t* clip; */
        uint16_t flag;
        int num_hits;
        int map_q;
        int error;
        int n_cigar;
        int n_alloc_cigar;
        int aln_len;
        int aln_len_alloc;
        uint32_t n_clip5;
        uint32_t n_clip3;
        /* uint64_t stop; */
        /* Per-record CIGAR counts (filled by parse_alignment, standard-plots.md 3.4):
         * m_bases = bases in M, = and X ops; *_ev = number of I / D operations;
         * *_bases = bases in them; s_bases = soft-clipped bases; x_ops = bases in X
         * ops. All zero for unmapped records and for a CIGAR that does not fit SEQ. */
        uint32_t m_bases;
        uint32_t i_ev;
        uint32_t i_bases;
        uint32_t d_ev;
        uint32_t d_bases;
        uint32_t s_bases;
        uint32_t x_ops;
        uint8_t reverse;
        uint8_t has_md_or_x;    /* 1 if the record has an MD tag or X ops (mismatches are known) */
};

/* Reader options. The batch limits and the -p cutoff are applied by the reader. */
struct reader_opts {
        int pool_threads;        /* decompression threads (0: none), one pool per file */
        uint64_t batch_records;  /* a chunk ends after this many records ... */
        uint64_t batch_bases;    /* ... or when it holds this many bases, whichever is first */
        uint64_t max_primary;    /* -p: stop after exactly this many primary records */
};

/* A batch of raw records as read from the file (WP4.1): the reader fills batches, the
 * workers convert them. Slots are allocated on first use and recycled, so a batch costs
 * nothing until records arrive; the data of a slot that held a very long record is
 * released by bam_batch_reset(). */
struct bam_batch {
        bam1_t **recs;
        int n;                   /* records in use */
        int alloc;               /* slots allocated (of the recs array) */
        int cap;                 /* most records one batch may hold */
};

/* One reader for every input: SAM, BAM, FASTA and FASTQ, plain or compressed, from a
 * file or from standard input ("-"). The type is taken from the content. */
struct sam_bam_file {
        samFile *in;
        bam_hdr_t* header;       /* NULL for an empty input */
        bam1_t *b;
        struct bam_batch *chunk; /* used by read_bam_chunk only */
        char* filename;
        hts_tpool *pool;
        struct htsThreadPool tp;
        uint64_t n_secondary;     /* skipped secondary alignments (0x100) */
        uint64_t n_supplementary; /* skipped supplementary alignments (0x800) */
        uint64_t n_primary;       /* primary records handed out so far */
        uint64_t batch_records;
        uint64_t batch_bases;
        uint64_t max_primary;
        uint8_t is_fastx;         /* FASTA or FASTQ: records carry no aln_data */
        uint8_t is_empty;         /* zero-byte or whitespace-only input */
        uint8_t has_sq;           /* the header has at least one @SQ line */
        uint8_t at_end;           /* nothing more to hand out */
        uint8_t partial;          /* -p cutoff reached and more primary records follow */
};


EXTERN int add_aln_data(struct tl_seq_buffer *sb);
EXTERN int clear_aln_data(struct tl_seq_buffer *sb);
EXTERN int free_aln_data(struct tl_seq_buffer *sb);

/* Batches: bam_batch_alloc(&b, cap) with cap = the record limit of one batch. */
EXTERN int bam_batch_alloc(struct bam_batch **batch, int cap);
EXTERN void bam_batch_reset(struct bam_batch *batch);
EXTERN void bam_batch_free(struct bam_batch *batch);
/* Reader side: read the next records into <batch> (it must be reset). Counts secondary
 * and supplementary records, skips empty FASTA/FASTQ records and applies the batch limits
 * and the -p cutoff. batch->n == 0 on return means that nothing more follows. Must be
 * called by one thread only. */
EXTERN int read_bam_batch(struct sam_bam_file *f_handle, struct bam_batch *batch);
/* Worker side: append the records of <batch> to <sb> (names, bases, qualities and, for
 * SAM/BAM, the alignment data). Touches no reader state, so any thread may call it. */
EXTERN int convert_bam_batch(const struct bam_batch *batch, int is_fastx, struct tl_seq_buffer *sb);
/* Read and convert in one step (single-threaded callers and tests). */
EXTERN int read_bam_chunk(struct sam_bam_file *f_handle, struct tl_seq_buffer *sb);
EXTERN int open_reads(struct sam_bam_file **file, char *filename, const struct reader_opts *opts);
EXTERN int open_bam(struct sam_bam_file **file, char *filename); /* defaults: no pool, no cutoff */
EXTERN int close_bam(struct sam_bam_file* f_handle);

#undef HTSGLUE_IMPORT
#undef EXTERN


#endif
