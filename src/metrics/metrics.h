#ifndef METRICS_H
#define METRICS_H

#include "seq/tld-seq.h"
#include <stdint.h>
#ifdef METRICS_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct seq_composition {
        uint32_t** data;
        int len;
        int L;
};


struct metrics {
        struct seq_composition* seq_comp;

        uint32_t min_len;
        uint32_t max_len;
};

EXTERN int metrics_alloc(struct metrics **metrics);
EXTERN void metrics_free(struct metrics *m);
EXTERN int get_metrics(struct tl_seq_buffer *sb, struct metrics *m);

EXTERN int debug_metrics_print(struct metrics *m);
#undef METRICS_IMPORT
#undef EXTERN


#endif
