#ifndef METRICS_H
#define METRICS_H

#include "seq/tld-seq.h"
#include <stdint.h>
#ifdef METRICS_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAPQUALBIN_gt30 0
#define MAPQUALBIN_lt30 1
#define MAPQUALBIN_lt10 2
#define MAPQUALBIN_UNMAP 3

struct mapqual_bins {
        uint8_t* map;
        char** description;
        int len;
        int n_bin;
};

struct error_composition {
        uint32_t** mis;
        uint32_t* ins;
        uint32_t* del;
        int len;
        int L;
        uint32_t n_mis;
        uint32_t n_ins;
        uint32_t n_del;
};

struct seq_composition {
        uint32_t** data;
        int L;
        uint32_t n_counts;
        int len;
};

struct qual_composition {
        uint32_t** data;
        int len;
        int L;
        uint32_t n_counts;
};

struct metrics {
        struct seq_composition** seq_comp_R1;
        struct qual_composition** qual_comp_R1;
        struct error_composition** error_comp_R1;


        struct seq_composition** seq_comp_R2;
        struct qual_composition** qual_comp_R2;
        struct error_composition** error_comp_R2;

        struct mapqual_bins* mapq_map;
        uint32_t min_len_R1;
        uint32_t max_len_R1;

        uint32_t min_len_R2;
        uint32_t max_len_R2;
        uint32_t report_max_len;
        uint8_t n_mapq_bins;
        uint8_t is_aligned;

        uint32_t n_R1_reads;
        uint32_t n_R2_reads;
        uint32_t n_paired;
        uint32_t n_proper_paired;
};

/* EXTERN int metrics_alloc(struct metrics **metrics); */
EXTERN int metrics_alloc(struct metrics **metrics, int report_max_len);
EXTERN void metrics_free(struct metrics *m);
EXTERN int get_metrics(struct tl_seq_buffer *sb, struct metrics *m);

/* EXTERN int debug_metrics_print(struct metrics *m); */
#undef METRICS_IMPORT
#undef EXTERN


#endif
