#ifndef HTSGLUE_H
#define HTSGLUE_H

#include "seq/tld-seq.h"
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
        int num_hits;
        int map_q;
        int error;
        uint64_t start;
        uint64_t stop;
};

struct sam_bam_file {
        samFile *in;
        bam_hdr_t* header;
        bam1_t *b;
        char* filename;

};


EXTERN int add_aln_data(struct tl_seq_buffer *sb);
EXTERN int clear_aln_data(struct tl_seq_buffer *sb);

EXTERN int read_bam_chunk(struct sam_bam_file *f_handle, struct tl_seq_buffer *sb);
EXTERN int open_bam(struct sam_bam_file **file, char *filename);
EXTERN int close_bam(struct sam_bam_file* f_handle);

#undef HTSGLUE_IMPORT
#undef EXTERN


#endif
