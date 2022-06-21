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

#define MAX_PST_MODEL_DEPTH 30

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
};

int pst_model_create(struct pst_model** model, struct tl_seq_buffer* sb);

#undef PST_IMPORT
#undef EXTERN

#endif
