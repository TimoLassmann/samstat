#ifndef HTSGLUE_H
#define HTSGLUE_H

#include "tld.h"
#include "sam.h"
#include "hts.h"
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

struct sam_bam_file {
        samFile *in;
        bam_hdr_t* header;
        bam1_t *b;
        char* filename;
        uint64_t n_secondary;     /* skipped secondary alignments (0x100) */
        uint64_t n_supplementary; /* skipped supplementary alignments (0x800) */
};


EXTERN int add_aln_data(struct tl_seq_buffer *sb);
EXTERN int clear_aln_data(struct tl_seq_buffer *sb);
EXTERN int free_aln_data(struct tl_seq_buffer *sb);

EXTERN int read_bam_chunk(struct sam_bam_file *f_handle, struct tl_seq_buffer *sb);
EXTERN int open_bam(struct sam_bam_file **file, char *filename);
EXTERN int close_bam(struct sam_bam_file* f_handle);

#undef HTSGLUE_IMPORT
#undef EXTERN


#endif
