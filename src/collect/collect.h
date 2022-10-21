#ifndef COLLECT_H
#define COLLECT_H

#include <stdint.h>

#ifdef COLLECT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif
struct tl_seq_buffer;
struct mapqual_bins;
struct plot_data;
struct stat_collection {
        /* one struct per plot...  */
        struct plot_data* base_comp_R1;
        struct plot_data* base_comp_R2;

        struct plot_data* qual_comp_R1;
        struct plot_data* qual_comp_R2;

        struct plot_data* mis_R1;
        struct plot_data* mis_R2;

        struct plot_data* ins_R1;
        struct plot_data* ins_R2;

        struct plot_data* del_R1;
        struct plot_data* del_R2;

        struct plot_data* len_dist_R1;
        struct plot_data* len_dist_R2;

        struct mapqual_bins* mapq_map;

        uint64_t n_read1;
        uint64_t n_read2;
        uint64_t n_paired;
        uint64_t n_proper_paired;
};


EXTERN int collect_stats(struct tl_seq_buffer *sb, struct stat_collection *s);
EXTERN int stat_collection_finalise(struct stat_collection *s);
EXTERN int stat_collection_alloc(struct stat_collection **stats);
EXTERN void stat_collection_free(struct stat_collection *s);
#undef COLLECT_IMPORT
#undef EXTERN


#endif
