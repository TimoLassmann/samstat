
#ifndef SAMSTAT_HEADER

#define SAMSTAT_HEADER 


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdlib.h>
#include <stdio.h>

#include <stdarg.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <getopt.h>
#include <string.h>
#include <math.h>

#include "tldevel.h"
#include "thr_pool.h"


struct parameters {
        char** infile; /**< @brief Names of input files. */
        char* outfile;
        int infiles;/**<  @brief Number of input files. */
        int quiet_flag;
        int num_query;/**< @brief Number of sequences to read at one time. */
        char* format;
        char* filter;
        char* train;
        char* exact5;
        char* messages;
        char* buffer;
        int gzipped;
        int bzipped;
        int dust;
        int sam;
        int fasta;
        int local_out;
};


struct parameters* interface(int argc, char *argv[]);
void free_param(struct parameters* param);


unsigned int nuc_code[256];
unsigned int rev_nuc_code[5];

#endif
