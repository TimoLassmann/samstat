#ifndef PARAM_H
#define PARAM_H

#include <stdint.h>
#ifdef PARAM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

/* process exit codes */
#define SAMSTAT_EXIT_FAILED 1 /* at least one input could not be processed */
#define SAMSTAT_EXIT_USAGE 2 /* bad command line */

/* Batch limits of the reader: a chunk is 4096 records or 8 Mbases, whichever comes first. */
#define SAMSTAT_BATCH_RECORDS 4096
#define SAMSTAT_BATCH_BASES 8000000

/* A record longer than this fails its input: it is not a read (a chromosome of a genome
 * FASTA, a corrupt record), and holding several of them at once would need many GB. The
 * longest sequencing reads are a few Mb. */
#define SAMSTAT_MAX_READ_LEN 50000000

/* The progress counter on a terminal is updated every this many reads. */
#define SAMSTAT_PROGRESS_STEP 1000000

/* name of the input that stands for standard input */
#define SAMSTAT_STDIN "-"

#define SAMSTAT_MAX_THREADS 64 /* -t above this is a usage error */

struct samstat_param {
        char** infile;
        char* outfile;
        char* outdir;
        uint8_t* no_output;     /* 1: report directory not writable, input is failed up front */
        int buffer_size;        /* reader batch: records */
        uint64_t batch_bases;   /* reader batch: bases */
        uint64_t max_read_len;  /* longer records fail the input */
        uint64_t progress_step; /* reads between progress counter updates */
        int n_infile;
        int n_threads;
        int verbose;
        uint8_t json;
        uint8_t multiqc;
        uint64_t top_n;
        int32_t report_max_len;
        double subsample;
        uint64_t seed;
        /* uint8_t pst; */
};

EXTERN int parse_param(int argc, char *argv[], struct samstat_param **param);
EXTERN  int param_init(struct samstat_param **param);
EXTERN void param_free(struct samstat_param *p);
/* Name that report files of input <id> are derived from: the input path as given, or
 * "stdin" for standard input ("-"). Valid as long as p is. */
EXTERN const char *samstat_output_stem(const struct samstat_param *p, int id);
/* Decompression threads for the file reader of one input (0: none). */
EXTERN int samstat_pool_threads(const struct samstat_param *p);
/* Worker threads that convert, parse and collect for one input file (0: everything runs
 * in the calling thread). With -t T: 0 for T = 1 and when there are at least T inputs
 * (one input per thread); 1 for T = 2; T - 1 - pool threads for T >= 3. */
EXTERN int samstat_worker_threads(const struct samstat_param *p);
/* Number of inputs processed at the same time: T if there are at least T inputs, else 1
 * (the inputs then run one after another, each with the whole thread budget). */
EXTERN int samstat_file_threads(const struct samstat_param *p);

#undef PARAM_IMPORT
#undef EXTERN


#endif
