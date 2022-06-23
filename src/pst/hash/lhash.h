#ifndef LHASH_H
#define LHASH_H

#ifdef LHASH_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#include <stdint.h>

EXTERN int init_lhash(void ** hash);

EXTERN int search_hash(void *hash,uint8_t* p, uint32_t len, uint32_t* counts);
EXTERN int insert_lhash(void *hash,uint8_t* p, uint32_t len);

EXTERN int delete_lhash(void* hash,uint8_t* p, uint32_t len);

EXTERN int reset_lhash(void* hash);
EXTERN void free_lhash(void* hash);

EXTERN void print_lhash(void* hash);

#undef LHASH_IMPORT
#undef EXTERN


#endif
