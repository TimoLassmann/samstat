#ifndef PARAM_H
#define PARAM_H

#ifdef PARAM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct samstat_param {
        char** infile; /**< @brief Names of input files. */
        char* outfile;
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
