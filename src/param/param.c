#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "htslib/hts_log.h"

#include "tldlite.h"
#define PARAM_IMPORT
#include "param.h"

#include "config.h"

#define OPT_VERBOSE 1
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
static int cap_threads_at_cpus(struct samstat_param *p);
static int same_file(const struct stat *a, const struct stat *b);
static int check_single_dash_long(int argc, char *argv[], const struct option *opts);

int print_help(char **argv);
static void print_usage(FILE *f, char **argv);

int parse_param(int argc, char *argv[], struct samstat_param **param)
{
        static const struct option long_options[] ={
                {"peek",required_argument,0,'p'},
                {"len",required_argument,0,'l'},   /* deprecated, ignored */
                {"dir",required_argument,0,'d'},
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
        struct samstat_param* p = NULL;
        int c;
        int help = 0;
        int version = 0;
        int in_options = 1;     /* an error here is a usage error: print the short usage */
        uint64_t v = 0;

        RUN(param_init(&p));
        RUN(check_single_dash_long(argc, argv, long_options));

        while (1){
                int option_index = 0;
                /* getopt_long, not getopt_long_only: long options need "--" (POSIX/GNU), so
                 * "-json" is "-j -s -o -n" and an error, not a silent --json */
                c = getopt_long (argc, argv,"p:l:d:t:hv",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case OPT_VERBOSE:
                        p->verbose = 1;
                        break;
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
                        /* owned copy: replaced by its realpath() later, freed in param_free() */
                        if(p->outdir){
                                MFREE(p->outdir);
                        }
                        p->outdir = strdup(optarg);
                        ASSERT(p->outdir != NULL, "out of memory");
                        break;
                case 'p':
                        RUN(parse_unsigned(optarg, "-p/--peek", UINT64_MAX, &v));
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
                        /* an option in long_options without a case above */
                        ERROR_MSG("internal error: unhandled option code %d", c);
                }
        }
        in_options = 0;
        /* htslib warns once per record for some problems (a reference name missing from the
         * header): gigabytes of log on a large file. Its errors are always shown, its
         * warnings only with --verbose. */
        hts_set_log_level(p->verbose ? HTS_LOG_WARNING : HTS_LOG_ERROR);

        if(help){
                print_help(argv);
                param_free(p);
                exit(EXIT_SUCCESS);
        }

        if(version){
                fprintf(stdout,"%s %s\n",PACKAGE_NAME,PACKAGE_VERSION);
                param_free(p);
                exit(EXIT_SUCCESS);
        }

        if( (argc- optind) <= 0){
                /* usage error: the short usage goes to stderr, the help (-h) to stdout */
                print_usage(stderr, argv);
                param_free(p);
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
        if(in_options){
                print_usage(stderr, argv);
        }
        param_free(p);
        return FAIL;
}

/* Long options need two dashes. Up to 2.3 (getopt_long_only) "-dir" and "-verbose" were long
 * options; read as short options they would silently mean "-d ir" and "-v" (print the version),
 * so such an argument is a usage error that names the two-dash form. "--" ends the options;
 * the argument of an option that takes one is skipped. */
static int check_single_dash_long(int argc, char *argv[], const struct option *opts)
{
        for(int i = 1; i < argc; i++){
                const char *a = argv[i];
                size_t n;

                if(strcmp(a, "--") == 0){
                        break;
                }
                if(a[0] != '-' || a[1] == 0){
                        continue;
                }
                if(a[1] == '-'){
                        /* --name VALUE: skip VALUE */
                        n = strcspn(a + 2, "=");
                        if(a[2 + n] == 0){
                                for(const struct option *o = opts; o->name; o++){
                                        if(o->has_arg == required_argument && strlen(o->name) == n && strncmp(o->name, a + 2, n) == 0){
                                                i++;
                                                break;
                                        }
                                }
                        }
                        continue;
                }
                if(a[2] == 0){
                        /* -p VALUE, -l VALUE, -d VALUE, -t VALUE: skip VALUE */
                        if(strchr("pldt", a[1]) != NULL){
                                i++;
                        }
                        continue;
                }
                n = strcspn(a + 1, "=");
                for(const struct option *o = opts; o->name; o++){
                        if(strlen(o->name) == n && strncmp(o->name, a + 1, n) == 0){
                                ERROR_MSG("Option %s: long options need two dashes (--%s)", a, o->name);
                        }
                }
        }
        return OK;
ERROR:
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
        RUN(cap_threads_at_cpus(p));
        if(p->top_n < 1){
                ERROR_MSG("Option -p/--peek: has to be >= 1.");
        }

        /* -d is checked before any input: a directory that cannot be created or written is a
         * usage error (exit 2), whether it exists already or not */
        if(p->outdir){
                char *ptr = NULL;
                if(make_dir_p(p->outdir) != OK){
                        ERROR_MSG("Option -d/--dir: cannot create output directory %s: %s", p->outdir, strerror(errno));
                }
                if(access(p->outdir, W_OK | X_OK) != 0){
                        ERROR_MSG("Option -d/--dir: cannot write to %s: %s; use -d DIR to write elsewhere.", p->outdir, strerror(errno));
                }
                ptr = realpath(p->outdir, NULL);
                if(!ptr){
                        ERROR_MSG("Could not obtain the full path for %s.", p->outdir);
                }
                MFREE(p->outdir);
                p->outdir = ptr;
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
        /* identical inputs (also under different names) are rejected in check_outputs() */
        RUN(check_outputs(p));

        return OK;
ERROR:
        return FAIL;
}

/* -t is the total thread budget: more threads than CPUs only add contention and memory (every
 * worker holds its own batches and counts), so a larger value is reduced to the number of
 * online CPUs. The output does not depend on -t. */
static int cap_threads_at_cpus(struct samstat_param *p)
{
        long n_cpu = sysconf(_SC_NPROCESSORS_ONLN);

        if(n_cpu < 1){
                n_cpu = 1;      /* unknown */
        }
        if((long) p->n_threads > n_cpu){
                if(p->verbose){
                        LOG_MSG("-t %d: only %ld CPUs online, using %ld threads.", p->n_threads, n_cpu, n_cpu);
                }
                p->n_threads = (int) n_cpu;
        }
        return OK;
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

/* The report files of one input: <stem>.samstat.html always, the others with --json and
 * --multiqc (the names stat_report.c, json_report.c and mqc_report.c write). */
static const char *const report_suffix[] = {
        ".samstat.html",
        ".samstat.json",
        ".samstat_general_mqc.json",
        ".samstat_mapq_mqc.json",
        ".samstat_length_mqc.json",
        ".samstat_basequal_mqc.json",
        ".samstat_error_mqc.json"
};
#define N_REPORT_SUFFIX ((int) (sizeof(report_suffix) / sizeof(report_suffix[0])))

static int report_written(const struct samstat_param *p, int k)
{
        if(k == 0){
                return 1;
        }
        if(k == 1){
                return p->json;
        }
        return p->multiqc;
}

/* the same file: also through a symbolic or hard link, or a different spelling of the path */
static int same_file(const struct stat *a, const struct stat *b)
{
        return a->st_dev == b->st_dev && a->st_ino == b->st_ino;
}

/* Up-front input and output checks, before any record is read.
 * Usage errors (exit 2), nothing is read:
 * - the same input given twice (also under another name of the same file);
 * - two inputs that map to the same report file; compared without regard to case, since
 *   the default file systems of macOS (APFS, HFS+) and Windows (NTFS) ignore it;
 * - a report file that is one of the inputs (it would be overwritten while being read).
 * Per-input failures (exit 1, the other inputs are processed), reported here once:
 * - an input that does not exist or cannot be read (said before anything about the
 *   output, so a missing input in a missing directory is not reported as "cannot write");
 * - without -d: the directory next to the input is not writable (-d DIR is checked in
 *   sanity_check_param()). */
static int check_outputs(struct samstat_param *p)
{
        char **key = NULL;
        struct stat *in_st = NULL;
        uint8_t *in_ok = NULL;
        char *path = NULL;
        int ret = FAIL;

        MMALLOC(key, sizeof(char*) * (size_t) p->n_infile);
        MMALLOC(in_st, sizeof(struct stat) * (size_t) p->n_infile);
        MMALLOC(in_ok, sizeof(uint8_t) * (size_t) p->n_infile);
        for(int i = 0; i < p->n_infile; i++){
                key[i] = NULL;
                in_ok[i] = 0;
        }
        /* the report location of every input: <directory, resolved>/<file name> */
        for(int i = 0; i < p->n_infile; i++){
                const char *in = samstat_output_stem(p, i);
                const char *slash = strrchr(in, '/');
                const char *base = slash ? slash + 1 : in;
                char *dir = NULL;
                char *rdir = NULL;
                size_t len;

                if(strcmp(p->infile[i], SAMSTAT_STDIN) != 0 && stat(p->infile[i], &in_st[i]) == 0){
                        in_ok[i] = 1;
                }
                if(p->outdir){
                        dir = strdup(p->outdir);
                }else if(slash){
                        len = (size_t) (slash - in);
                        dir = strndup(in, len ? len : 1);   /* "/x" -> "/" */
                }else{
                        dir = strdup(".");
                }
                ASSERT(dir != NULL, "out of memory");
                rdir = realpath(dir, NULL);
                len = strlen(rdir ? rdir : dir) + strlen(base) + 2;
                key[i] = malloc(len);
                if(key[i] == NULL){
                        free(dir);
                        free(rdir);
                        ERROR_MSG("out of memory");
                }
                snprintf(key[i], len, "%s/%s", rdir ? rdir : dir, base);
                free(dir);
                free(rdir);
        }
        /* usage errors */
        for(int i = 0; i < p->n_infile; i++){
                for(int j = i + 1; j < p->n_infile; j++){
                        int dup;
                        if(in_ok[i] && in_ok[j]){
                                dup = same_file(&in_st[i], &in_st[j]);
                        }else{
                                dup = strcmp(p->infile[i], p->infile[j]) == 0;
                        }
                        if(dup){
                                ERROR_MSG("Duplicated input files: %s %s", p->infile[i], p->infile[j]);
                        }
                }
        }
        for(int i = 0; i < p->n_infile; i++){
                for(int j = i + 1; j < p->n_infile; j++){
                        if(strcasecmp(key[i], key[j]) == 0){
                                ERROR_MSG("Output collision: %s and %s would both write %s.samstat.html (file names are compared without regard to case); use -d DIR or rename an input.",
                                          p->infile[i], p->infile[j], key[i]);
                        }
                }
        }
        for(int i = 0; i < p->n_infile; i++){
                for(int k = 0; k < N_REPORT_SUFFIX; k++){
                        struct stat st;
                        size_t len;

                        if(!report_written(p, k)){
                                continue;
                        }
                        len = strlen(key[i]) + strlen(report_suffix[k]) + 1;
                        path = malloc(len);
                        ASSERT(path != NULL, "out of memory");
                        snprintf(path, len, "%s%s", key[i], report_suffix[k]);
                        if(stat(path, &st) == 0){
                                for(int j = 0; j < p->n_infile; j++){
                                        if(in_ok[j] && same_file(&st, &in_st[j])){
                                                ERROR_MSG("Output collision: the report %s of %s would overwrite the input %s; use -d DIR or rename an input.",
                                                          path, p->infile[i], p->infile[j]);
                                        }
                                }
                        }
                        free(path);
                        path = NULL;
                }
        }
        /* per-input failures */
        for(int i = 0; i < p->n_infile; i++){
                const char *in = samstat_output_stem(p, i);
                const char *slash = strrchr(in, '/');
                const char *base = slash ? slash + 1 : in;
                char *dir = NULL;
                size_t len;

                if(strcmp(p->infile[i], SAMSTAT_STDIN) != 0 && access(p->infile[i], R_OK) != 0){
                        if(errno == ENOENT){
                                fprintf(stderr, "samstat: error: input %s does not exist.\n", p->infile[i]);
                        }else{
                                fprintf(stderr, "samstat: error: cannot read input %s: %s.\n", p->infile[i], strerror(errno));
                        }
                        p->no_output[i] = 1;
                        continue;
                }
                if(p->outdir){
                        continue;
                }
                if(slash){
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
                free(dir);
        }
        ret = OK;
ERROR:
        free(path);
        if(key){
                for(int i = 0; i < p->n_infile; i++){
                        free(key[i]);
                }
                MFREE(key);
        }
        if(in_st){
                MFREE(in_st);
        }
        if(in_ok){
                MFREE(in_ok);
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
        fprintf(stdout,"  -t, --threads N    total threads, 1-64, at most the number of CPUs [4]\n");
        fprintf(stdout,"      --json         also write <input>.samstat.json\n");
        fprintf(stdout,"      --multiqc      also write MultiQC files <input>.samstat_*_mqc.json\n");
        fprintf(stdout,"      --verbose      debug output\n");
        fprintf(stdout,"  -h, --help / -v, --version\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Input type (SAM, BAM, CRAM, FASTA or FASTQ, plain or gzip) is detected from the\n");
        fprintf(stdout,"content. CRAM needs its reference, found through REF_PATH or REF_CACHE (no\n");
        fprintf(stdout,"network). A file name of - reads standard input (once); its reports are\n");
        fprintf(stdout,"named stdin.samstat.*.\n");
        fprintf(stdout,"Exit codes: 0 ok, 1 one or more inputs failed, 2 usage error (bad option,\n");
        fprintf(stdout,"-d DIR not writable, colliding report names; nothing is read).\n");
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


/* After a usage error: on stderr, so that it does not mix with output a caller reads. */
static void print_usage(FILE *f, char **argv)
{
        char* tmp = NULL;
        const char* name = PACKAGE_NAME;

        if(tld_get_filename(argv[0], &tmp) == OK){
                name = tmp;
        }
        fprintf(f, "Usage: %s [options] <file|-> [file ...]\n", name);
        fprintf(f, "Try '%s --help' for more information.\n", name);
        if(tmp){
                MFREE(tmp);
        }
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
        p->outdir = NULL;
        p->no_output = NULL;
        p->top_n = UINT64_MAX;
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
                if(p->outdir){
                        MFREE(p->outdir);
                }

                MFREE(p);
        }

}
