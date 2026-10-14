#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

#include "misc/misc.h"
#include "tld.h"
#define PARAM_IMPORT
#include "param.h"

#include "config.h"

#define OPT_VERBOSE 1
/* #define OPT_SUBSAMPLE 2 */
#define OPT_SEED 3
#define OPT_PLOTEND 4
#define OPT_PLOTNORM 5
#define OPT_JSON 6



int sanity_check_param(struct samstat_param *p);
static int parse_unsigned(const char *arg, const char *opt, uint64_t max, uint64_t *out);
int detect_file_type(char *filename, int* type);

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
                        {"len",required_argument,0,'l'},
                        {"dir",required_argument,0,'d'},
                        /* {"sub",required_argument,0,OPT_SUBSAMPLE}, */
                        {"seed",required_argument,0,OPT_SEED},
                        {"plotend",0,0,OPT_PLOTEND},
                        {"norm",0,0,OPT_PLOTNORM},
                        {"json",0,0,OPT_JSON},
                        {"nthreads",required_argument,0,'t'},
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
                case OPT_PLOTEND:
                        p->collect_end = 1;
                        break;

                case OPT_PLOTNORM:
                        p->norm_len_plot = 1;
                        break;
                case OPT_JSON:
                        p->json = 1;
                        break;
                case 't':
                        RUN(parse_unsigned(optarg, "-t", INT32_MAX, &v));
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
                        RUN(parse_unsigned(optarg, "-l/-len", INT32_MAX, &v));
                        p->report_max_len = (int32_t) v;
                        /* param->local_out = 1; */
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

        galloc(&p->file_type, p->n_infile);
        for(int i = 0; i < p->n_infile;i++){
                p->file_type[i] = FILE_TYPE_UNKNOWN;
        }
        for(int i = 0 ; i < p->n_infile;i++){
                int t = -1;
                RUN(detect_file_type(p->infile[i], &t));
                p->file_type[i] = t;
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
                if(tld_dir_exists(p->outdir) == FAIL){
                        ERROR_MSG("Option -d/-dir: output directory %s does not exist",p->outdir);
                }
                char *ptr = NULL;
                ptr = realpath(p->outdir, NULL);
                if(!ptr){
                        ERROR_MSG("Could not obtain the full path for %s.", p->outdir);
                }
                p->outdir = ptr;

                /* int l = strnlen(p->outdir, 1024); */
                /* if(p->outdir[l-1] == '/'){ */
                /*         p->outdir[l-1] = 0; */
                /* } */
                /* char *symlinkpath = "/tmp/symlink/file"; */
                /* char actualpath [PATH_MAX+1]; */
                /* ptr = realpath(symlinkpath, actualpath); */
        }
        if(p->subsample < 0.0){
                ERROR_MSG("Option -sub: has to be greater than 0.0");
        }else if(p->subsample > 1.0){
                ERROR_MSG("Option -sub: has to be smaller than 1.0");
        }

        /* check report length */
        if(p->report_max_len <= 0){
                ERROR_MSG("Option -l/-len: the report length has to be > 0");
        }else if(p->report_max_len > 10000){
                ERROR_MSG("Option -l/-len: the report length is > 10000, the maximum currently supported.");
        }else if(p->report_max_len > 500){
                WARNING_MSG("Option -l/-len: the report length is > 500. This will slow down samstat considerably.");
        }
        /* check input files  */

        int n_unknown = 0;
        int n_sam = 0;
        int n_fasta = 0;
        for(int i = 0; i < p->n_infile;i++){
                if(p->file_type[i] == FILE_TYPE_UNKNOWN){
                        n_unknown++;
                }else if(p->file_type[i] == FILE_TYPE_FASTAQ){
                        n_fasta++;
                }else if(p->file_type[i] == FILE_TYPE_SAMBAM){
                        n_sam++;
                }else {
                        ERROR_MSG("File type %d not recognized - this should never happen", (int) p->file_type[i]);
                }
        }

        if(n_unknown){
                for(int i = 0; i < p->n_infile;i++){
                        if(p->file_type[i] == FILE_TYPE_UNKNOWN){
                                LOG_MSG("   unrecognized file name extension: %s", p->infile[i]);
                        }
                }
                ERROR_MSG("Could not recognize the type of %d input file(s) by their extension (supported: .sam .bam .fa .fasta .fna .fq .fastq, and .gz variants of the fasta/fastq ones).", n_unknown);
        }
        if(n_fasta && n_sam){
                WARNING_MSG("Found both fasta/fastq and sam/bam files.");
        }
        /* check for identical files */
        for(int i = 0; i < p->n_infile;i++){
                for(int j = i+1; j < p->n_infile;j++){
                        if(strncmp(p->infile[i], p->infile[j], 512) == 0){
                                ERROR_MSG("Duplicated input files: %s %s",p->infile[i], p->infile[j]);
                        }
                }
        }


        return OK;
ERROR:
        return FAIL;
}

static int has_suffix_nocase(const char *name, const char *suffix)
{
        size_t n = strlen(name);
        size_t m = strlen(suffix);
        if(n <= m){
                return 0;
        }
        return strcasecmp(name + n - m, suffix) == 0;
}

int detect_file_type(char *filename, int* type)
{
        static const char* sambam_suffix[] = {".sam", ".bam", NULL};
        static const char* fastaq_suffix[] = {
                ".fa", ".fasta", ".fna", ".fq", ".fastq",
                ".fa.gz", ".fasta.gz", ".fna.gz", ".fq.gz", ".fastq.gz",
                NULL};

        *type = FILE_TYPE_UNKNOWN;
        for(int i = 0; sambam_suffix[i]; i++){
                if(has_suffix_nocase(filename, sambam_suffix[i])){
                        *type = FILE_TYPE_SAMBAM;
                        return OK;
                }
        }
        for(int i = 0; fastaq_suffix[i]; i++){
                if(has_suffix_nocase(filename, fastaq_suffix[i])){
                        *type = FILE_TYPE_FASTAQ;
                        return OK;
                }
        }
        /* a missing file is reported when the file is processed (exit code 1) */
        return OK;
}

int print_help(char **argv )
{
        const char usage[] = " <file1> <file2> ....";
        char* tmp = NULL;
        RUN(tld_get_filename(argv[0], &tmp));


        fprintf(stdout,"\n");
        fprintf(stdout,"%s (v%s)\n", tmp, PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2011,2022 Timo Lassmann\n");
        fprintf(stdout,"\n");


        fprintf(stdout,"Usage: %s [-options] %s\n\n",tmp,usage);

        fprintf(stdout,"Options:\n\n");
        /* fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-p","Report over-represented patterns." ,"[]"  ); */
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-d/-dir","Output directory.","[]");
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",22-3,"","NOTE: by default SAMStat will place reports in the same directory as the input files.","");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-p/-peek","Report stats on exactly the first <n> sequences (n >= 1).","[unlimited]");
        /* fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-l/-len","Report stats on the first <n> nucleotides." ,"[250]"  ); */
        /* fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"--sub","Sub-sample of reads." ,"[1.0]"  ); */
        /* fprintf(stdout,"%*s%-*s  %s %s\n",3,"",22-3,"","e.g. \"--sub 0.2\" would report stats on a random 20% selection of reads." ,""  ); */
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-t","Number of threads." ,"[4]"  );
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",22-3,"","will only be used when multiple input files are present." ,""  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"--plotend","Add base and quality plots relative to the read ends." ,"[]"  );
         fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"--norm","Add plots summarising statistics in 1% read length intervals." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"--json","Also write machine-readable statistics to <input>.samstat.json." ,"[]"  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"--seed","Random number seed." ,"[0]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"--verbose","Enables verbose output." ,"[]"  );

        fprintf(stdout,"\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-h/-help","Prints help message." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-v/-version","Prints version information." ,"[]"  );

        fprintf(stdout,"\n");

        fprintf(stdout,"Examples:\n\n");

        fprintf(stdout,"   %s ~/tmp/test.bam \n",tmp);
        fprintf(stdout,"   - will write a QC report to: ~/tmp/test.bam.samstat.html\n\n");

        fprintf(stdout,"   %s ~/tmp/test.bam -p 1000000 -d ~/samstat_reports/\n",tmp);
        fprintf(stdout,"   - will write a QC report based on the top 1000000 sequences to: ~/samstat_reports/test.bam.samstat.html\n\n")
;

        fprintf(stdout,"   %s ~/tmp/*.bam -p 1000000 -d ~/samstat_reports/\n",tmp);
        fprintf(stdout,"   - will write QC reports for all BAM files in ~/tmp, based on the top 1000000 sequences to: ~/samstat_reports/test.bam.samstat.html\n\n")
;

        /* fprintf(stdout,"\n"); */
        /* fprintf(stdout,"\n"); */
        /* fprintf(stdout,"Please cite:\n"); */

        /* fprintf(stdout,"  Lassmann, Timo, Yoshihide Hayashizaki, and Carsten O. Daub.\n"); */

        /* fprintf(stdout,"  \"SAMStat: monitoring biases in next generation sequencing data.\"\n"); */
        /* fprintf(stdout,"  Bioinformatics 27.1 (2011): 130-131. \n"); */
        /* fprintf(stdout,"  https://doi.org/10.1093/bioinformatics/btq614\n"); */

        /* fprintf(stdout,"\n"); */

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
        p->buffer_size = 10000;
        p->verbose = 0;
        p->n_infile = 0;
        p->n_threads = 4;
        p->collect_end = 0;
        p->norm_len_plot = 0;
        p->json = 0;
        p->infile = NULL;
        p->outfile = NULL;
        p->outdir = NULL;
        p->file_type = NULL;
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
                if(p->file_type){
                        gfree(p->file_type);
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
