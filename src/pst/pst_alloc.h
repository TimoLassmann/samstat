#ifndef PST_ALLOC_H
#define PST_ALLOC_H

#include <stdint.h>

#ifdef PST_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct pst;
struct pst_node;
struct fpst;

EXTERN int init_pst(struct pst** pst, double min_error, double expected_error, int len, int L);
EXTERN int reset_pst(struct pst* p);
EXTERN void free_pst(struct pst *p);

EXTERN int alloc_fpst(struct fpst** fast_pst, int size, int L);
EXTERN int resize_fpst(struct fpst* f, int L);
EXTERN void free_fpst(struct fpst* f);

EXTERN void free_pst_node(struct pst_node* n, int L);
EXTERN int alloc_node(struct pst_node** node,uint8_t* string,int len, int L);

#undef PST_ALLOC_IMPORT
#undef EXTERN


#endif
