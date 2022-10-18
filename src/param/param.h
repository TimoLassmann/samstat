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
        char* outdir;
        uint8_t* file_type;
        int buffer_size;
        int n_infile;
        int verbose;
        uint64_t top_n;
        int32_t report_max_len;
        double subsample;
        uint64_t seed;
        /* uint8_t pst; */
};

EXTERN int parse_param(int argc, char *argv[], struct samstat_param **param);
EXTERN  int param_init(struct samstat_param **param);
EXTERN void param_free(struct samstat_param *p);

#undef PARAM_IMPORT
#undef EXTERN


#endif
