#ifndef TOOLS_H
#define TOOLS_H

#ifdef TOOLS_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#include <stdint.h>
struct rng_state;

EXTERN int sample_wo_replacement(int pop_size, int sample_size,struct rng_state* rng, int* sel);
EXTERN int get_file_size(char *filename, uint64_t *size);



#undef TOOLS_IMPORT
#undef EXTERN


#endif
