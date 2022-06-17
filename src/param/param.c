#include <getopt.h>
#include <unistd.h>

#include "core/tld-core.h"
#include "tld.h"
#define PARAM_IMPORT
#include "param.h"

#include "../config.h"

void usage();

int parse_param(int argc, char *argv[], struct samstat_param **param)
{
        struct samstat_param* p = NULL;
        int c;
        int help = 0;
        int version = 0;

        param_init(&p);


        while (1){
                static struct option long_options[] ={
                        {"help",0,0,'h'},
                        {"version",0,0,'v'},
                        {"log",required_argument,0,'l'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"hvl",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case 'h':
                        help = 1;
                        break;
                case 'v':
                        version = 1;
                        break;
                case 'l':
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
                usage();
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
                usage();
                param_free(p);
                /* free_param(param); */
                exit(EXIT_SUCCESS);
        }
        LOG_MSG("%d %d", argc, optind);


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
        *param = p;
        return OK;
ERROR:
        param_free(p);
        return FAIL;
}

void usage()
{
        fprintf(stderr, "\n%s %s, Copyright (C) 2022 Timo Lassmann\n",PACKAGE_NAME, PACKAGE_VERSION);
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   samstat <file1> <file2>  ...  \n\n");
        fprintf(stderr, "SAMstat will produce a summary file (html) for each input file named\n <original filename>.samstat.html.\n");
        fprintf(stderr, "\n");
}


int param_init(struct samstat_param **param)
{
        struct samstat_param* p = NULL;

        MMALLOC(p, sizeof(struct samstat_param));
        p->buffer_size = 10000;
        p->quiet = 0;
        p->n_infile = 0;
        p->infile = NULL;
        p->outfile = NULL;
        p->file_type = NULL;

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
