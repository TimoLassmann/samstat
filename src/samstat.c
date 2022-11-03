#include "tld.h"
#include "sambamparse/sam_bam_parse.h"
#include "htslib/sam.h"
#include "htsinterface/htsglue.h"
#include "param/param.h"
#include "collect/collect.h"
#include "report/stat_report.h"
#include "tools/tools.h"
#include <stdint.h>
#include <inttypes.h>

#if HAVE_PTHREADS
#include "thread/thread_data.h"
#endif

/* int detect_file_type(char *filename, int* type); */
int process_sam_bam_file(struct samstat_param *p, int id);
int process_fasta_fastq_file(struct samstat_param *p, int id);

int debug_seq_buffer_print(struct tl_seq_buffer *sb);
void *samstat_worker(void *data);

int main(int argc, char *argv[])
{
        struct samstat_param* param = NULL;

        RUN(parse_param(argc, argv, &param));

#if HAVE_PTHREADS
        struct thread_data* td = NULL;
        RUN(thread_data_init(&td, param));

        for(long t = 0; t < td->n_threads; t++){
                int rc = pthread_create(&td->threads[t], NULL, samstat_worker, (void *)td);
                if(rc){
                        ERROR_MSG("ERROR; return code from pthread_create() is %d\n", rc);
                }
        }
        for(long t = 0; t < td->n_threads; t++) {
                void *status;
                int rc = pthread_join(td->threads[t], &status);
                if(rc){
                        ERROR_MSG("ERROR; return code from pthread_join() is %d\n", rc);
                }
                if(status == PTHREAD_CANCELED){
                        fprintf(stdout,"Main: completed join with thread %ld having a status of %ld (cancelled)\n",t,(long)status);
                }
        }

        thread_data_free(td);
#else
        for(int i = 0 ; i < param->n_infile;i++){
                if(param->verbose){
                        LOG_MSG("Processing: %s", param->infile[i]);
                }
                int t = -1;
                t = param->file_type[i];
                if(t == FILE_TYPE_FASTAQ){
                        RUN(process_fasta_fastq_file(param,i));
                }else if(t == FILE_TYPE_SAMBAM){
                        RUN(process_sam_bam_file(param,i));
                }
        }
#endif
        param_free(param);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

#if HAVE_PTHREADS
void *samstat_worker(void *data)
{
        struct thread_data* td= (struct thread_data*) data;

        struct samstat_param* param = td->p;

        int run = 1;

        while(run){
                int target_id = -1;
                /* find work  */
                 if(pthread_mutex_lock(&td->lock) != 0){
                          ERROR_MSG("Can't get lock");
                 }

                 for(int i = 0; i < param->n_infile;i++){
                         if(td->active[i] == 0){
                                 td->active[i] = 1;
                                 target_id = i ;
                                 break;
                         }
                 }

                 if(param->verbose && target_id != -1){
                         LOG_MSG("Working on %s",param->infile[target_id]);
                 }

                 if(pthread_mutex_unlock(&td->lock) != 0){
                         ERROR_MSG("Can't get lock");
                 }

                 if(target_id == -1){
                         run = 0;
                 }else{
                         int t = -1;
                         t = param->file_type[target_id];
                         if(t == FILE_TYPE_FASTAQ){
                                 RUN(process_fasta_fastq_file(param,target_id));
                         }else if(t == FILE_TYPE_SAMBAM){
                                 RUN(process_sam_bam_file(param,target_id));
                         }
                 }
        }

        pthread_exit(0);
        return NULL;
ERROR:
        LOG_MSG("Thread %d done due to error.");
        pthread_exit(0);
        return NULL;
}
#endif

int process_sam_bam_file(struct samstat_param* p, int id)
{
        struct sam_bam_file* f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct alphabet* a = NULL;
        /* struct metrics* metrics = NULL; */

        struct stat_collection* s = NULL;
        uint64_t n_read = 0;
        uint64_t old_n_read = 0;

        ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);

        RUN(stat_collection_alloc(&s));
        RUN(stat_collection_config_additional_plot(s,p));
        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));
        add_aln_data(sb);
        RUN(create_alphabet(&a, NULL,TLALPHABET_NOAMBIGUOUS_DNA));

        sb->data = a;

        s->is_aligned = 1;

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

                RUN(collect_stats(sb, s));

                n_read += sb->num_seq;
                if(p->verbose){
                        if(n_read >= old_n_read + 1000000){
                                LOG_MSG("Processed %"PRId64" sequences", n_read);
                                old_n_read = n_read;
                        }
                }
                if(n_read > p->top_n && sb->num_seq == sb->malloc_num){
                        if(p->verbose){
                                LOG_MSG("Stopping because more than %"PRId64" sequences were read.", n_read);
                        }
                        s->is_partial_report = 1;
                        break;
                }
                clear_aln_data(sb);
                reset_tl_seq_buffer(sb);
        }
        RUN(close_bam(f_handle));

        RUN(stat_report(s,p, id));
        stat_collection_free(s);

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
        struct stat_collection* s = NULL;
        uint64_t n_read = 0;

        ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);

        RUN(open_fasta_fastq_file(&f_handle, p->infile[id], TLSEQIO_READ));
        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));
        RUN(stat_collection_alloc(&s));
        RUN(stat_collection_config_additional_plot(s,p));

        while(1){
                RUN(read_fasta_fastq_file(f_handle, &sb, p->buffer_size));
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

                if(sb->num_seq == 0){
                        break;
                }

                RUN(collect_stats(sb, s));

                n_read += sb->num_seq;
                if(p->verbose){
                        LOG_MSG("Processed %"PRId64" sequences", n_read);
                }
                if(n_read > p->top_n){
                        if(p->verbose){
                                LOG_MSG("Stopping because more than %"PRId64" sequences were read.", n_read);
                        }
                        s->is_partial_report = 1;
                        break;
                }
                reset_tl_seq_buffer(sb);
        }

        RUN(stat_report(s,p, id));
        stat_collection_free(s);
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
