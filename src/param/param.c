#include <getopt.h>
#include <unistd.h>

#include "core/tld-core.h"
#include "tld.h"
#define PARAM_IMPORT
#include "param.h"

#include "../config.h"

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
                        {"len",required_argument,0,'l'},
                        {"help",0,0,'h'},
                        {"version",0,0,'v'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"hvpl:",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case 'p':
                        p->pst = 1;
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
                //fprintf(stdout,"%s %s\n",PACKAGE_NAME,PACKAGE_VERSION);
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
        param_free(p);
        return FAIL;
}

int sanity_check_param(struct samstat_param *p)
{

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


int print_help(char **argv )
{
        const char usage[] = " <file1> <file2> ....";
        char* tmp = NULL;

        RUN(tlfilename(argv[0], &tmp));
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",tmp,usage);

        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-p","Report over-represented patterns." ,"[]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",22-3,"-l/-len","Report stats on the first <l> nucleotides." ,"[500]"  );


        fprintf(stdout,"\n");
        return OK;
ERROR:
        return FAIL;
}


int param_init(struct samstat_param **param)
{
        struct samstat_param* p = NULL;

        MMALLOC(p, sizeof(struct samstat_param));
        p->buffer_size = 100000;
        p->quiet = 0;
        p->n_infile = 0;
        p->infile = NULL;
        p->outfile = NULL;
        p->file_type = NULL;
        p->pst = 0;
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
                MFREE(p);
        }

}
