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
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>

#if HAVE_PTHREADS
#include "thread/thread_data.h"
#endif

/* Read input <id> (SAM, BAM, FASTA or FASTQ, from a file or standard input) and report. */
int process_reads_file(struct samstat_param *p, int id);

int debug_seq_buffer_print(struct tl_seq_buffer *sb);
void *samstat_worker(void *data);

/* Run the analysis for input file <id> depending on its type. */
static int process_file(struct samstat_param *p, int id);

int main(int argc, char *argv[])
{
        struct samstat_param* param = NULL;
        int n_failed = 0;
        int ret = EXIT_SUCCESS;

        if(parse_param(argc, argv, &param) != OK){
                return SAMSTAT_EXIT_USAGE;
        }

#if HAVE_PTHREADS
        struct thread_data* td = NULL;
        RUN(thread_data_init(&td, param));

        int n_created = 0;
        for(long t = 0; t < td->n_threads; t++){
                int rc = pthread_create(&td->threads[n_created], NULL, samstat_worker, (void *)td);
                if(rc){
                        WARNING_MSG("pthread_create() failed with return code %d.", rc);
                }else{
                        n_created++;
                }
        }
        if(n_created == 0){
                /* no worker thread could be started - do the work here */
                samstat_worker((void *) td);
        }
        for(int t = 0; t < n_created; t++) {
                void *status;
                int rc = pthread_join(td->threads[t], &status);
                if(rc){
                        WARNING_MSG("pthread_join() failed with return code %d.", rc);
                }
        }
        for(int i = 0; i < param->n_infile; i++){
                if(!td->ok[i]){
                        n_failed++;
                }
        }
        thread_data_free(td);
#else
        for(int i = 0 ; i < param->n_infile;i++){
                if(param->verbose){
                        LOG_MSG("Processing: %s", param->infile[i]);
                }
                if(process_file(param, i) != OK){
                        WARNING_MSG("Failed to process file: %s", param->infile[i]);
                        n_failed++;
                }
        }
#endif
        param_free(param);
        if(n_failed){
                WARNING_MSG("%d of the input files could not be processed.", n_failed);
                ret = SAMSTAT_EXIT_FAILED;
        }
        return ret;
#if HAVE_PTHREADS
ERROR:
        param_free(param);
        return SAMSTAT_EXIT_FAILED;
#endif
}

static int process_file(struct samstat_param *p, int id)
{
        if(p->no_output[id]){
                /* already reported up front by param.c: fail this input silently */
                goto ERROR;
        }
        RUN(process_reads_file(p, id));
        return OK;
ERROR:
        return FAIL;
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
                 }else if(process_file(param, target_id) == OK){
                         td->ok[target_id] = 1;
                 }else{
                         /* an error in one file must not stop the remaining files */
                         WARNING_MSG("Failed to process file: %s", param->infile[target_id]);
                 }
        }
        return NULL;
ERROR:
        LOG_MSG("Worker thread done due to error.");
        return NULL;
}
#endif

int process_reads_file(struct samstat_param* p, int id)
{
        struct sam_bam_file* f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct alphabet* a = NULL;
        struct stat_collection* s = NULL;
        struct reader_opts opts;
        uint64_t n_read = 0;
        uint64_t old_n_read = 0;
        int is_stdin = (strcmp(p->infile[id], SAMSTAT_STDIN) == 0);
        int ret = FAIL;

        if(!is_stdin){
                ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);
        }

        opts.pool_threads = samstat_pool_threads(p);
        opts.batch_records = (uint64_t) p->buffer_size;
        opts.batch_bases = p->batch_bases;
        opts.max_primary = p->top_n;

        RUN(open_reads(&f_handle, p->infile[id], &opts));

        RUN(stat_collection_alloc(&s));
        RUN(stat_collection_config_additional_plot(s,p));
        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));
        if(!f_handle->is_fastx){
                add_aln_data(sb);
        }
        RUN(create_alphabet(&a, NULL,TLALPHABET_NOAMBIGUOUS_DNA));

        sb->data = a;

        /* NOTE: replaced by has_sq / mapped reads in the collector's finalise (WP3.1a) */
        s->is_aligned = f_handle->is_fastx ? 0 : 1;

        while(1){
                RUN(read_bam_chunk(f_handle, sb));
                sb->L = TL_SEQ_BUFFER_DNA;
                sb->base_quality_offset = 33; /* Can I get this from sam - don't need to part of the standard */

                if(sb->num_seq == 0){
                        break;
                }

                if(!f_handle->is_fastx){
                        for(int i = 0; i < sb->num_seq;i++){
                                parse_alignment(sb->sequences[i]);
                                fix_orientation(sb->sequences[i]);
                        }
                }

                RUN(collect_stats(sb, s));

                n_read += sb->num_seq;
                if(p->verbose){
                        if(n_read >= old_n_read + 1000000){
                                LOG_MSG("Processed %"PRId64" sequences", n_read);
                                old_n_read = n_read;
                        }
                }
                if(!f_handle->is_fastx){
                        clear_aln_data(sb);
                }
                reset_tl_seq_buffer(sb);
        }
        s->n_secondary = f_handle->n_secondary;
        s->n_supplementary = f_handle->n_supplementary;
        s->is_partial_report = f_handle->partial;
        ret = close_bam(f_handle);
        f_handle = NULL; /* closed - must not be closed again in the cleanup below */
        RUN(ret);
        ret = FAIL;

        if(n_read == 0){
                WARNING_MSG("%s: the input is empty (no reads), the report will be empty.", p->infile[id]);
        }
        if(s->is_partial_report){
                WARNING_MSG("%s: partial report, based on the first %"PRIu64" reads only.", p->infile[id], n_read);
        }
        RUN(stat_report(s,p, id));
        ret = OK;
ERROR:
        /* shared by the success and the error path */
        if(f_handle){
                close_bam(f_handle);
        }
        stat_collection_free(s);
        if(sb){
                if(sb->data){
                        free_alphabet((struct alphabet*) sb->data);
                        sb->data = NULL;
                }
                free_aln_data(sb);
                free_tl_seq_buffer(sb);
        }
        return ret;
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
