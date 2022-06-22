#ifndef PST_H
#define PST_H

#include "seq/tld-seq.h"
#include "tld.h"
#include <stdint.h>
#ifdef PST_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAX_PST_MODEL_DEPTH 12

struct pst;
/* #include "pst_structs.h" */
struct pst_model{
        void* h;                /* count hash  */
        double* counts_l;       /* how many kmers present */
        int L;                  /* Alphabet size  */
        uint32_t depth;              /* depth of model - how many previous states are examined
                                     to determine the p of next symbol */
        struct pst* pst;        /* The actual pst model  */
        double min_error;
        double gamma;
        uint32_t n_seq;
};


int pst_model_alloc(struct pst_model **model);
int pst_model_add_seq(struct pst_model *m, struct tl_seq_buffer *sb);
int pst_model_create(struct pst_model *m);

#undef PST_IMPORT
#undef EXTERN

#endif
