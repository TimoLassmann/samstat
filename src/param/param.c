#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "tldlite.h"
#define PARAM_IMPORT
#include "param.h"

#include "config.h"

#define OPT_VERBOSE 1
/* #define OPT_SUBSAMPLE 2 */
#define OPT_SEED 3
#define OPT_PLOTEND 4
#define OPT_PLOTNORM 5
#define OPT_JSON 6
#define OPT_MULTIQC 7



int sanity_check_param(struct samstat_param *p);
static int make_dir_p(const char *path);
static int check_outputs(struct samstat_param *p);
static void deprecated_option(const char *opt);
static int parse_unsigned(const char *arg, const char *opt, uint64_t max, uint64_t *out);
static void batch_limits_from_env(struct samstat_param *p);

int print_help(char **argv);

int parse_param(int argc, char *argv[], struct samstat_param **param)
{
        struct samstat_param* p = NULL;
        int c;
        int help = 0;
        int version = 0;
        uint64_t v = 0;

        RUN(param_init(&p));

        while (1){
                static struct option long_options[] ={
                        {"peek",required_argument,0,'p'},
                        {"len",required_argument,0,'l'},   /* deprecated, ignored */
                        {"dir",required_argument,0,'d'},
                        /* {"sub",required_argument,0,OPT_SUBSAMPLE}, */
                        {"seed",required_argument,0,OPT_SEED},
                        {"plotend",0,0,OPT_PLOTEND},
                        {"norm",0,0,OPT_PLOTNORM},
                        {"json",0,0,OPT_JSON},
                        {"multiqc",0,0,OPT_MULTIQC},
                        {"threads",required_argument,0,'t'},
                        {"nthreads",required_argument,0,'t'}, /* alias */
                        {"verbose",0,0,OPT_VERBOSE},
                        {"help",0,0,'h'},
                        {"version",0,0,'v'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"p:l:d:t:hv",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case OPT_VERBOSE:
                        p->verbose = 1;
                        break;
                /* case OPT_SUBSAMPLE: */
                /*         p->subsample = atof(optarg); */
                /*         break; */
                case OPT_SEED:
                        deprecated_option("--seed");
                        break;
                case OPT_PLOTEND:
                        deprecated_option("--plotend");
                        break;
                case OPT_PLOTNORM:
                        deprecated_option("--norm");
                        break;
                case OPT_JSON:
                        p->json = 1;
                        break;
                case OPT_MULTIQC:
                        p->multiqc = 1;
                        break;
                case 't':
                        RUN(parse_unsigned(optarg, "-t/--threads", INT32_MAX, &v));
                        ASSERT(v <= SAMSTAT_MAX_THREADS, "Option -t/--threads: at most %d threads (got %s)", SAMSTAT_MAX_THREADS, optarg);
                        p->n_threads = (int) v;
                        break;
                case 'd':
                        p->outdir = optarg;
                        break;
                case 'p':
                        RUN(parse_unsigned(optarg, "-p/-peek", UINT64_MAX, &v));
                        p->top_n = v;
                        break;
                case 'h':
                        help = 1;
                        break;
                case 'v':
                        version = 1;
                        break;
                case 'l':
                        deprecated_option("-l");
                        break;
                case '?':
                        /* getopt has already printed the offending option */
                        goto ERROR;
                        break;
                default:
                        fprintf(stderr,"default\n\n\n\n");
                        abort ();
                }
        }
        if(help){
                print_help(argv);
                param_free(p);
                /* free_param(param); */
                exit(EXIT_SUCCESS);
        }

        if(version){
                fprintf(stdout,"%s %s\n",PACKAGE_NAME,PACKAGE_VERSION);
                param_free(p);
                exit(EXIT_SUCCESS);
        }

        if( (argc- optind) <= 0){
                print_help(argv);
                param_free(p);
                /* free_param(param); */
                exit(SAMSTAT_EXIT_USAGE);
        }

        MMALLOC(p->infile,sizeof(char*)* (argc-optind));

        c = 0;
        while (optind < argc){
                p->infile[c] =  argv[optind++];
                c++;
        }
        p->n_infile = c;

        galloc(&p->no_output, p->n_infile);
        for(int i = 0; i < p->n_infile;i++){
                p->no_output[i] = 0;
        }

        RUN(sanity_check_param(p));
        *param = p;
        return OK;
ERROR:
        WARNING_MSG("Something went wrong with parsing the command line options.");

        /* WARNING_MSG("Something went wrong with parsing the command line options."); */
        /* print_help(argv); */
        param_free(p);
        return FAIL;
}

/* Strict unsigned integer parsing: digits only (no sign, no blanks, no
 * trailing characters), no overflow, value <= max. */
static int parse_unsigned(const char *arg, const char *opt, uint64_t max, uint64_t *out)
{
        char *end = NULL;
        unsigned long long v;

        ASSERT(arg != NULL, "Option %s: missing value", opt);
        ASSERT(arg[0] >= '0' && arg[0] <= '9', "Option %s: '%s' is not a non-negative integer", opt, arg);
        errno = 0;
        v = strtoull(arg, &end, 10);
        ASSERT(*end == 0, "Option %s: '%s' is not a non-negative integer", opt, arg);
        ASSERT(errno != ERANGE && v <= max, "Option %s: value '%s' is too large", opt, arg);
        *out = (uint64_t) v;
        return OK;
ERROR:
        return FAIL;
}

int sanity_check_param(struct samstat_param *p)
{

        if(p->n_threads < 1){
                ERROR_MSG("Option -t (num threads) : has to be > 0");
        }
        if(p->top_n < 1){
                ERROR_MSG("Option -p/-peek: has to be >= 1.");
        }



        if(p->outdir){
                char *ptr = NULL;
                if(make_dir_p(p->outdir) != OK){
                        ERROR_MSG("Option -d/-dir: cannot create output directory %s: %s", p->outdir, strerror(errno));
                }
                ptr = realpath(p->outdir, NULL);
                if(!ptr){
                        ERROR_MSG("Could not obtain the full path for %s.", p->outdir);
                }
                p->outdir = ptr;
        }
        if(p->subsample < 0.0){
                ERROR_MSG("Option -sub: has to be greater than 0.0");
        }else if(p->subsample > 1.0){
                ERROR_MSG("Option -sub: has to be smaller than 1.0");
        }

        /* check input files: the type is taken from the content when the file is opened,
         * so an unrecognised input fails only that input */
        {
                int n_stdin = 0;
                for(int i = 0; i < p->n_infile;i++){
                        if(strcmp(p->infile[i], SAMSTAT_STDIN) == 0){
                                n_stdin++;
                        }
                }
                if(n_stdin > 1){
                        ERROR_MSG("Standard input (-) can be given only once.");
                }
        }
        /* check for identical files */
        for(int i = 0; i < p->n_infile;i++){
                for(int j = i+1; j < p->n_infile;j++){
                        if(strncmp(p->infile[i], p->infile[j], 512) == 0){
                                ERROR_MSG("Duplicated input files: %s %s",p->infile[i], p->infile[j]);
                        }
                }
        }

        RUN(check_outputs(p));

        return OK;
ERROR:
        return FAIL;
}

static void deprecated_option(const char *opt)
{
        fprintf(stderr, "samstat: warning: %s is deprecated and ignored\n", opt);
}

/* mkdir -p */
static int make_dir_p(const char *path)
{
        char *tmp = NULL;
        size_t n = strlen(path);
        struct stat st;

        if(n == 0){
                errno = ENOENT;
                return FAIL;
        }
        tmp = strdup(path);
        if(!tmp){
                return FAIL;
        }
        while(n > 1 && tmp[n-1] == '/'){
                tmp[--n] = 0;
        }
        for(char *q = tmp + 1; *q; q++){
                if(*q == '/'){
                        *q = 0;
                        if(mkdir(tmp, 0777) != 0 && errno != EEXIST){
                                int e = errno;
                                free(tmp);
                                errno = e;
                                return FAIL;
                        }
                        *q = '/';
                }
        }
        if(mkdir(tmp, 0777) != 0 && errno != EEXIST){
                int e = errno;
                free(tmp);
                errno = e;
                return FAIL;
        }
        free(tmp);
        if(stat(path, &st) != 0){
                return FAIL;
        }
        if(!S_ISDIR(st.st_mode)){
                errno = ENOTDIR;
                return FAIL;
        }
        return OK;
}

/* Up-front output checks, before any record is read:
 * - the directory a report goes to must be writable: otherwise that input is
 *   marked as failed (exit 1);
 * - two inputs must not map to the same report file (usage error, exit 2). */
static int check_outputs(struct samstat_param *p)
{
        char **key = NULL;
        int ret = FAIL;

        MMALLOC(key, sizeof(char*) * (size_t) p->n_infile);
        for(int i = 0; i < p->n_infile; i++){
                key[i] = NULL;
        }
        for(int i = 0; i < p->n_infile; i++){
                const char *in = samstat_output_stem(p, i);
                const char *slash = strrchr(in, '/');
                const char *base = slash ? slash + 1 : in;
                char *dir = NULL;
                char *rdir = NULL;
                size_t len;

                if(p->outdir){
                        dir = strdup(p->outdir);
                }else if(slash){
                        len = (size_t) (slash - in);
                        dir = strndup(in, len ? len : 1);   /* "/x" -> "/" */
                }else{
                        dir = strdup(".");
                }
                ASSERT(dir != NULL, "out of memory");
                if(access(dir, W_OK | X_OK) != 0){
                        fprintf(stderr, "samstat: error: cannot write %s/%s.samstat.html: %s; use -d DIR to write elsewhere.\n",
                                dir, base, strerror(errno));
                        p->no_output[i] = 1;
                }
                rdir = realpath(dir, NULL);
                len = strlen(rdir ? rdir : dir) + strlen(base) + 2;
                key[i] = malloc(len);
                ASSERT(key[i] != NULL, "out of memory");
                snprintf(key[i], len, "%s/%s", rdir ? rdir : dir, base);
                free(dir);
                free(rdir);
        }
        for(int i = 0; i < p->n_infile; i++){
                for(int j = i + 1; j < p->n_infile; j++){
                        if(strcmp(key[i], key[j]) == 0){
                                ERROR_MSG("Output collision: %s and %s would both write %s.samstat.html; use -d DIR or rename an input.",
                                          p->infile[i], p->infile[j], key[i]);
                        }
                }
        }
        ret = OK;
ERROR:
        if(key){
                for(int i = 0; i < p->n_infile; i++){
                        free(key[i]);
                }
                MFREE(key);
        }
        return ret;
}

/* The only place the test-only limits are read. SAMSTAT_BATCH_RECORDS, SAMSTAT_BATCH_BASES,
 * SAMSTAT_MAX_READ_LEN and SAMSTAT_PROGRESS_STEP can LOWER the built-in value (to cut a small
 * test input into many batches, or to test a limit without a huge input). They can never raise
 * it: memory is allocated in proportion to these values, so a large value must not get past
 * here. They are not a user option and are not documented; anything but a positive integer is
 * ignored. */
static uint64_t env_lower(const char *name, uint64_t builtin)
{
        const char *e = getenv(name);
        char *end = NULL;
        unsigned long long v;

        if(!e || e[0] < '0' || e[0] > '9'){
                return builtin;
        }
        errno = 0;
        v = strtoull(e, &end, 10);
        if(*end != 0 || errno == ERANGE || v == 0 || v >= builtin){
                return builtin;
        }
        return (uint64_t) v;
}

static void batch_limits_from_env(struct samstat_param *p)
{
        p->buffer_size = (int) env_lower("SAMSTAT_BATCH_RECORDS", SAMSTAT_BATCH_RECORDS);
        p->batch_bases = env_lower("SAMSTAT_BATCH_BASES", SAMSTAT_BATCH_BASES);
        p->max_read_len = env_lower("SAMSTAT_MAX_READ_LEN", SAMSTAT_MAX_READ_LEN);
        p->progress_step = env_lower("SAMSTAT_PROGRESS_STEP", SAMSTAT_PROGRESS_STEP);
}

const char *samstat_output_stem(const struct samstat_param *p, int id)
{
        if(strcmp(p->infile[id], SAMSTAT_STDIN) == 0){
                return "stdin";
        }
        return p->infile[id];
}

int samstat_pool_threads(const struct samstat_param *p)
{
        /* -t is the total thread budget: 1 or 2 threads leave nothing for decompression;
         * files that fill all threads on their own run without a pool */
        if(p->n_threads < 3 || p->n_infile >= p->n_threads){
                return 0;
        }
        return MACRO_MAX(1, p->n_threads / 4);
}

int samstat_worker_threads(const struct samstat_param *p)
{
        /* The calling thread reads. T >= n_files: one file per thread, each inline.
         * Otherwise the files run one after another and spend all T threads: the reader,
         * P decompression threads (none below T = 3) and the rest as workers. */
        if(p->n_threads < 2 || p->n_infile >= p->n_threads){
                return 0;
        }
        return p->n_threads - 1 - samstat_pool_threads(p);
}

int samstat_file_threads(const struct samstat_param *p)
{
        /* threads that each work on a different file at the same time */
        if(p->n_infile >= p->n_threads){
                return p->n_threads;
        }
        return 1;
}

int print_help(char **argv )
{
        char* tmp = NULL;
        RUN(tld_get_filename(argv[0], &tmp));

        fprintf(stdout,"%s %s - quality control for SAM/BAM/CRAM and FASTA/FASTQ files\n", PACKAGE_NAME, PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Usage: %s [options] <file|-> [file ...]\n\n", tmp);
        fprintf(stdout,"  -d, --dir DIR      write reports to DIR (created if missing) [next to input]\n");
        fprintf(stdout,"  -p, --peek N       analyse only the first N reads (not random) [all]\n");
        fprintf(stdout,"  -t, --threads N    total threads, 1-64 [4]\n");
        fprintf(stdout,"      --json         also write <input>.samstat.json\n");
        fprintf(stdout,"      --multiqc      also write MultiQC files <input>.samstat_*_mqc.json\n");
        fprintf(stdout,"      --verbose      debug output\n");
        fprintf(stdout,"  -h, --help / -v, --version\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Input type (SAM, BAM, CRAM, FASTA or FASTQ, plain or gzip) is detected from the\n");
        fprintf(stdout,"content. CRAM needs its reference, found through REF_PATH or REF_CACHE (no\n");
        fprintf(stdout,"network). A file name of - reads standard input (once); its reports are\n");
        fprintf(stdout,"named stdin.samstat.*.\n");
        fprintf(stdout,"Exit codes: 0 ok, 1 one or more inputs failed, 2 usage error.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Examples:\n");
        fprintf(stdout,"  %s sample.bam\n", tmp);
        fprintf(stdout,"  %s -d qc/ run/*.bam run/*.fastq.gz\n", tmp);
        fprintf(stdout,"\n");

        MFREE(tmp);
        return OK;
ERROR:
        if(tmp){
                MFREE(tmp);
        }
        return FAIL;
}


int param_init(struct samstat_param **param)
{
        struct samstat_param* p = NULL;

        MMALLOC(p, sizeof(struct samstat_param));
        batch_limits_from_env(p);
        p->verbose = 0;
        p->n_infile = 0;
        p->n_threads = 4;
        p->json = 0;
        p->multiqc = 0;
        p->infile = NULL;
        p->outfile = NULL;
        p->outdir = NULL;
        p->no_output = NULL;
        p->top_n = UINT64_MAX;
        p->subsample = 1.0;
        /* p->pst = 0; */
        p->report_max_len = 250;
        p->seed = 0;
        *param = p;
        return OK;
ERROR:
        param_free(p);
        return FAIL;
}


void param_free(struct samstat_param *p)
{

        if(p){
                if(p->no_output){
                        gfree(p->no_output);
                }
                if(p->infile){
                        MFREE(p->infile);
                }
                /* if(p->outdir){ */
                        /* MFREE(p->outdir); */
                /* } */

                MFREE(p);
        }

}
