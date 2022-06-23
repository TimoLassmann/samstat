#ifndef MUR_H
#define MUR_H

#include <stdint.h>
#include <stdlib.h>

#ifdef MUR_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed);

#undef MUR_IMPORT
#undef EXTERN


#endif
