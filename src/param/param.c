#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "misc/misc.h"
#include "tld.h"
#define PARAM_IMPORT
#include "param.h"

#include "../config.h"

#define OPT_VERBOSE 1

int sanity_check_param(struct samstat_param *p);
int detect_file_type(char *filename, int* type);

int print_help(char **argv);

int parse_param(int argc, char *argv[], struct samstat_param **param)
{
        struct samstat_param* p = NULL;
        int c;
        int help = 0;
        int version = 0;

        param_init(&p);

        while (1){
                static struct option long_options[] ={
                        {"peek",required_argument,0,'p'},
                        {"len",required_argument,0,'l'},
                        {"dir",required_argument,0,'d'},
                        {"verbose",0,0,OPT_VERBOSE},
                        {"help",0,0,'h'},
                        {"version",0,0,'v'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"p:l:d:hv",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case OPT_VERBOSE:
                        p->verbose = 1;
                        break;
                case 'd':
                        p->outdir = optarg;
                        break;
                case 'p':
                        p->top_n = atoi(optarg);
                        break;
                case 'h':
                        help = 1;
                        break;
                case 'v':
                        version = 1;
                        break;
                case 'l':
                        p->report_max_len = atoi(optarg);
                        /* param->local_out = 1; */
                        break;
                case '?':
                        exit(1);
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
                exit(EXIT_SUCCESS);
        }
        /* LOG_MSG("%d %d", argc, optind); */


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

int sanity_check_param(struct samstat_param *p)
{
        if(p->top_n < 10){
                ERROR_MSG("Option -p/-peek: has to be >= 10.");
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
        /* check report length */
        if(p->report_max_len <= 0){
                ERROR_MSG("Option -l/-len: the report length has to be > 0");
        }else if(p->report_max_len > 500){
                WARNING_MSG("Option -l/-len: the report length is > 500. This will slow down samstat considerably.");
        }else if(p->report_max_len > 10000){
                ERROR_MSG("Option -l/-len: the report length is > 10000, the maximum currently supported.");
        }
        /* check input files  */

        int n_unknown = 0;
        int n_sam = 0;
        int n_fasta = 0;
        for(int i = 0; i < p->n_infile;i++){
                if(p->file_type[i] == FILE_TYPE_UNKNOWN){
                        n_unknown++;
                }else if(p->file_type[i] == FILE_TYPE_FASTA){
                        n_fasta++;
                }else if(p->file_type[i] == FILE_TYPE_SAMBAM){
                        n_sam++;
                }else {
                        ERROR_MSG("File type %d not recognized - this should never happen");
                }
        }

        if(p->n_infile == n_unknown){
                LOG_MSG("Could not recognize the type (fasta / sam / bam) of any input file by their file name extension:");
                for(int i = 0; i < p->n_infile;i++){
                        LOG_MSG("   %s", p->infile[i]);
                }

                ERROR_MSG("Quitting");


        }else{
                if(n_fasta && n_sam){
                        WARNING_MSG("Found both fasta and sam/bam files.");
                        for(int i = 0; i < p->n_infile;i++){
                                if(p->file_type[i] != FILE_TYPE_UNKNOWN){
                                        LOG_MSG("   %s", p->infile[i]);
                                }
                        }
                }
        }


        return OK;
ERROR:
        return FAIL;
}

int detect_file_type(char *filename, int* type)
{
        if(tld_file_exists(filename)== OK){
                tld_str f = tld_char_to_str(filename);
                tld_str sam_suffix = tld_char_to_str(".sam");
                tld_str bam_suffix = tld_char_to_str(".bam");
                tld_str fasta_suffix = tld_char_to_str(".fasta");
                tld_str fasta_suffix_short = tld_char_to_str(".fa");
                tld_str fastq_suffix = tld_char_to_str(".fastq");
                tld_str fastq_suffix_2 = tld_char_to_str(".fq");

                tld_str fasta_suffix_gz = tld_char_to_str(".fasta.gz");
                tld_str fasta_suffix_gz2 = tld_char_to_str(".fna.gz");
                tld_str fastq_suffix_gz = tld_char_to_str(".fastq.gz");

                if(tld_suffix_match(f, sam_suffix)){
                        *type = FILE_TYPE_SAMBAM;
                }else if(tld_suffix_match(f, bam_suffix)){
                        *type = FILE_TYPE_SAMBAM;
                }else if(tld_suffix_match(f, fasta_suffix)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fasta_suffix_short)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fastq_suffix)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fastq_suffix_2)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fasta_suffix_gz)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fasta_suffix_gz2)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fastq_suffix_gz)){
                        *type = FILE_TYPE_FASTAQ;
                }else{
                        *type = FILE_TYPE_UNKNOWN;
                }

        }else{
                *type = FILE_TYPE_UNKNOWN;
                ERROR_MSG("File %s not found!!!!", filename);
        }

        return OK;
ERROR:
        return FAIL;
}

/* int print_kalign_header(void) */
/* { */
/*         fprintf(stdout,"\n"); */
/*         fprintf(stdout,"SAMStat (%s)\n", PACKAGE_VERSION); */
/*         fprintf(stdout,"\n"); */
/*         fprintf(stdout,"Copyright (C) 2022 Timo Lassmann\n"); */
/*         fprintf(stdout,"\n"); */
/*         fprintf(stdout,"Please cite:\n"); */

/*         /\*        fprintf(stdout,"  Kalign 3: multiple sequence alignment of large data sets */
/* Timo Lassmann */
/* Bioinformatics, btz795, https://doi.org/10.1093/bioinformatics/btz795 */
/*         *\/ */
/*         fprintf(stdout,"  Lassmann, Timo.\n"); */
/*         fprintf(stdout,"  \"Kalign 3: multiple sequence alignment of large data sets.\"\n"); */
/*         fprintf(stdout,"  Bioinformatics (2019) \n"); */
/*         fprintf(stdout,"  https://doi.org/10.1093/bioinformatics/btz795\n"); */
/*         fprintf(stdout,"\n"); */


/*         return OK; */
/* } */


int print_help(char **argv )
{
        const char usage[] = " <file1> <file2> ....";
        char* tmp = NULL;

        RUN(tlfilename(argv[0], &tmp));


        fprintf(stdout,"\n");
        fprintf(stdout,"%s (%s)\n", tmp, PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2011,2022 Timo Lassmann\n");
        fprintf(stdout,"\n");



        fprintf(stdout,"Usage: %s [-options] %s\n\n",tmp,usage);

        fprintf(stdout,"Options:\n\n");
        /* fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-p","Report over-represented patterns." ,"[]"  ); */
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-d/-dir","Output directory.","[]");
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",22-3,"","NOTE: by default SAMStat will place reports in the same directory as the input files.","");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-p/-peek","Report stats only on the first <n> sequences.","[unlimited]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-l/-len","Report stats on the first <n> nucleotides." ,"[500]"  );
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
        p->buffer_size = 1000000;
        p->verbose = 0;
        p->n_infile = 0;
        p->infile = NULL;
        p->outfile = NULL;
        p->outdir = NULL;
        p->file_type = NULL;
        p->top_n = UINT64_MAX;
        /* p->pst = 0; */
        p->report_max_len = 500;
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
                if(p->outdir){
                        MFREE(p->outdir);
                }

                MFREE(p);
        }

}
