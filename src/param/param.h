#ifndef PARAM_H
#define PARAM_H

#include <stdint.h>
#ifdef PARAM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define FILE_TYPE_SAMBAM 0
#define FILE_TYPE_FASTAQ 1
#define FILE_TYPE_UNKNOWN 2

struct samstat_param {
        char** infile;
        char* outfile;
        uint8_t* file_type;
        int buffer_size;
        int n_infile;
        int quiet;
};

EXTERN int parse_param(int argc, char *argv[], struct samstat_param **param);
EXTERN  int param_init(struct samstat_param **param);
EXTERN void param_free(struct samstat_param *p);

#undef PARAM_IMPORT
#undef EXTERN


#endif