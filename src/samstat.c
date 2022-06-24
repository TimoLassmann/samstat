
#include "core/tld-core.h"
#include "seq/tld-seq.h"
#include "tld.h"
#include "sambamparse/sam_bam_parse.h"
#include "htslib/sam.h"
#include "htsinterface/htsglue.h"
#include "param/param.h"
#include "metrics/metrics.h"
#include "report/report.h"
#include "pst.h"
#include<stdint.h>

    /* int detect_file_type(char *filename, int* type); */
int process_sam_bam_file(struct samstat_param *p, int id);
int process_fasta_fastq_file(struct samstat_param *p, int id);

int debug_seq_buffer_print(struct tl_seq_buffer *sb);

int main(int argc, char *argv[])
{
        struct samstat_param* param = NULL;

        parse_param(argc, argv, &param);

        /* fprintf(stdout,"Hello world\n"); */
        /* fprintf(stdout,"HTS version:%d\n", HTS_VERSION); */
        int* g = NULL;

        galloc(&g,100);
        for(int i = 0; i < 100;i++){
                g[i] = i;
        }
        gfree(g);

        for(int i = 0 ; i < param->n_infile;i++){
                int t = -1;
                /* RUN(detect_file_type(param->infile[i], &t)); */
                t = param->file_type[i];
                if(t == FILE_TYPE_FASTAQ){
                        process_fasta_fastq_file(param,i);
                }else if(t == FILE_TYPE_SAMBAM){
                        process_sam_bam_file(param,i);
                }
        }

        param_free(param);
        /* exit(0); */
        /* struct sam_bam_file* f_handle = NULL; */
/*         struct tl_seq_buffer* sb = NULL; */

/*         RUN(alloc_tl_seq_buffer(&sb, 20)); */
/*         add_aln_data(sb); */

/*         RUN(open_bam(&f_handle, argv[1])); */
/*         while(1){ */
/*                 RUN(read_bam_chunk(f_handle, sb)); */
/*                 if(sb->num_seq == 0){ */
/*                         break; */
/*                 } */
/*                 for(int i = 0 ; i < sb->num_seq;i++){ */
/*                         struct aln_data* a = (struct aln_data*) sb->sequences[i]->data; */
/*                         fprintf(stdout,"%s\n%s\n", TLD_STR( sb->sequences[i]->name), */
/*                                 TLD_STR()b->sequences[i]->seq */
/*                                 /\* TLD_STR(a->md), *\/ */
/*                                 /\* a->error, *\/ */
/*                                 /\* a->reverse *\/ */
/*                                 ); */

/*                         parse_alignment(sb->sequences[i]); */
/*                         /\* get_readable_cigar(a, NULL); *\/ */
/*                         char* cig = NULL; */
/*                         get_readable_cigar(a, &cig); */
/*                         fprintf(stdout,"Cigar: %s\n",cig); */
/*                         fprintf(stdout,"   MD: %s\n",TLD_STR(a->md)); */
/*                         gfree(cig); */

/*                         fprintf(stdout, "%s (genome)\n",a->genome); */
/*                         for(int i = 0; i < a->aln_len;i++){ */
/*                                 if(a->genome[i] == a->read[i]){ */
/*                                         fprintf(stdout,"|"); */
/*                                 }else{ */
/*                                         fprintf(stdout," "); */
/*                                 } */
/*                         } */
/*                         fprintf(stdout,"\n"); */

/*                         fprintf(stdout, "%s (read)\n",a->read); */

/*                 } */
/*                 break; */
/* /\* sb->num_seq = 0; *\/ */
/*                 reset_tl_seq_buffer(sb); */
/*         } */
        /* RUN(close_bam(f_handle)); */


        /* clear_aln_data(sb); */
        /* free_tl_seq_buffer(sb); */
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int process_sam_bam_file(struct samstat_param* p, int id)
{
        struct sam_bam_file* f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct alphabet* a = NULL;
        /* struct pst_model* m = NULL; */
        struct metrics* metrics = NULL;
        /* p->buffer_size = 1000; */

        ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);

        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));
        add_aln_data(sb);
        RUN(create_alphabet(&a, NULL,TLALPHABET_NOAMBIGUOUS_DNA));

        sb->data = a;

        /* if(p->pst){ */
        /*         pst_model_alloc(&m); */
        /* } */


        RUN(metrics_alloc(&metrics, p->report_max_len));
        metrics->is_aligned = 1;
        RUN(open_bam(&f_handle, p->infile[id]));
        while(1){

                RUN(read_bam_chunk(f_handle, sb));
                sb->L = TL_SEQ_BUFFER_DNA;
                sb->base_quality_offset = 33; /* Can I get this from sam - don't need to part of the standard */

                if(sb->num_seq == 0){
                        break;
                }

                for(int i = 0; i < sb->num_seq;i++){
                        parse_alignment(sb->sequences[i]);
                        fix_orientation(sb->sequences[i]);
                }
                RUN(get_metrics(sb,metrics));

                /* if(p->pst && m->n_seq < 1000000){ */
                /*         pst_model_add_seq(m, sb); */
                /* } */

                /* LOG_MSG("L: %d",sb->L); */
                //debug_seq_buffier_print(sb);
                /* break; */
                /* struct pst_model* m = NULL; */
                /* pst_model_create(&m, sb); */
                LOG_MSG("Processed %d sequences", sb->num_seq);
                /* exit(0); */
                clear_aln_data(sb);
                reset_tl_seq_buffer(sb);
                /* exit(0); */
        }
        RUN(close_bam(f_handle));

        /* if(p->pst){ */
        /*         pst_model_create(m); */
        /* } */

        /* RUN(debug_metrics_print(metrics)); */
        create_report(metrics, p,id);
        metrics_free(metrics);
        if(sb->data){
                a = sb->data;
                free_alphabet(a);
        }
        free_aln_data(sb);
        free_tl_seq_buffer(sb);
        return OK;
ERROR:
        return FAIL;
}

int process_fasta_fastq_file(struct samstat_param* p, int id)
{
        struct file_handler *f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct alphabet* a = NULL;
        /* struct pst_model* m = NULL; */
        struct metrics* metrics = NULL;

        ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);

        RUN(open_fasta_fastq_file(&f_handle, p->infile[id], TLSEQIO_READ));
        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));

        RUN(metrics_alloc(&metrics, p->report_max_len));
        metrics->is_aligned = 0;

        /* if(p->pst){ */
        /*         pst_model_alloc(&m); */
        /* } */

        while(1){

                RUN(read_fasta_fastq_file(f_handle, &sb, p->buffer_size));

                /* LOG_MSG("%d", sb->base_quality_offset); */
                /* detect_format(sb); */
                if(!sb->data){
                        LOG_MSG("Setting alphabet");
                        if(sb->L == TL_SEQ_BUFFER_DNA){
                                RUN(create_alphabet(&a, NULL,TLALPHABET_NOAMBIGUOUS_DNA ));
                        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                                RUN(create_alphabet(&a, NULL,TLALPHABET_NOAMBIGIOUS_PROTEIN));
                        }else{
                                ERROR_MSG("Bio polymer size %d not recognized", sb->L);
                        }
                        sb->data = a;
                }
                LOG_MSG("%d BQ offset ", sb->base_quality_offset);
                //total_r+= sb->num_seq;
                LOG_MSG("Finished reading chunk: found %d ",sb->num_seq);


                if(sb->num_seq == 0){
                        break;
                }
                RUN(get_metrics(sb,metrics));

                /* LOG_MSG("starting pst"); */
                /* if(p->pst && m->n_seq < 1000000){ */
                /*         pst_model_add_seq(m, sb); */
                /* } */
                /* pst_model_create(&m, sb); */

                /* LOG_MSG("L: %d",sb->L); */
                //debug_seq_buffer_print(sb);
                /* for(int i = 0; i < sb->num_seq;i++){ */
                /*         fprintf(stdout, "%s\n", sb->sequences[i]->seq); */
                /* } */
                reset_tl_seq_buffer(sb);
        }
        /* if(p->pst){ */
        /*         pst_model_create(m); */
        /* } */

        /* RUN(debug_metrics_print(metrics)); */
        create_report(metrics, p,id);
        metrics_free(metrics);
        if(sb->data){
                a = sb->data;
                free_alphabet(a);
        }
        free_tl_seq_buffer(sb);
        RUN(close_seq_file(&f_handle));
        return OK;
ERROR:
        return  FAIL;
}

int debug_seq_buffer_print(struct tl_seq_buffer *sb)
{
        for(int i = 0 ; i < sb->num_seq;i++){
                struct aln_data* a = (struct aln_data*) sb->sequences[i]->data;
                fprintf(stdout,"%s\n%s\n", TLD_STR( sb->sequences[i]->name),
                        TLD_STR(sb->sequences[i]->seq)
                        /* TLD_STR(a->md), */
                        /* a->error, */
                        /* a->reverse */
                        );
                if(a){
                        parse_alignment(sb->sequences[i]);
                        /* get_readable_cigar(a, NULL); */
                        char* cig = NULL;
                        get_readable_cigar(a, &cig);
                        fprintf(stdout,"Cigar: %s\n",cig);
                        fprintf(stdout,"   MD: %s\n",TLD_STR(a->md));
                        gfree(cig);

                        fprintf(stdout, "%s (genome)\n",a->genome);
                        for(int i = 0; i < a->aln_len;i++){
                                if(a->genome[i] == a->read[i]){
                                        fprintf(stdout,"|");
                                }else{
                                        fprintf(stdout," ");
                                }
                        }
                        fprintf(stdout,"\n");

                        fprintf(stdout, "%s (read)\n",a->read);
                }
        }
        return OK;
}

/*
  # include <string.h>
  # include <stdio.h>
  int main(int argc, char *argv[]) {
  int k;
  if (argc == 1) return 1;
  k = str2bwt(argv[1]);
  printf("%.*s$%s\n", k, argv[1], &argv[1][k]);
  return 0;
  }
  int str2bwt(char *s) {
  int i,c=0,k=0,l=strlen(s);
  for(i=l-1;i>=0;--i) {
  int j,r=0,a=s[i];
  memmove(&s[i],&s[i+1],k);
  s[i+k]=a;
  for(j=i;j<i+k;++j) r+=(s[j]<=a);
  for(;j<l;++j) r+=(s[j]<a);
  k=r+1,c=a;
  }
  return k; // sentinel
  }
*/
