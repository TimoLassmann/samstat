#ifndef PST_STRUCTS_H
#define PST_STRUCTS_H
#include <stdint.h>
struct pst_node{
        struct pst_node** next;
        double* probability;
        double* value;
        uint32_t* label;
        int label_len;
};

struct fpst{
        double** prob;
        double** value;
        int** links;
        void* data;
        int l;                  /* length */
        int m;                  /* what is malloced */
};


struct pst {
        struct fpst* fpst_root;
        struct pst_node* root;
        double* background;
        int max_observed_len;
        double p_min;
        double gamma_min;
        double epsilon;
        int L;
        int len;
        int num_nodes;
};


#endif
