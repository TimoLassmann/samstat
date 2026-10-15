#include "tldlite.h"
#include "sambamparse/sam_bam_parse.h"
#include "htslib/sam.h"
#include "htsinterface/htsglue.h"
#include "param/param.h"
#include "collect/collect.h"
#include "report/stat_report.h"
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <stdio.h>

#if HAVE_PTHREADS
#include "thread/thread_data.h"
#endif

_Static_assert(SAMSTAT_MAX_READ_LEN == READER_MAX_READ_LEN, "one read length limit");

/* Read input <id> (SAM, BAM, FASTA or FASTQ, from a file or standard input) and report. */
int process_reads_file(struct samstat_param *p, int id);

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

/* ----------------------------------------------------------------------------------
 * One input file (WP4.1). The thread that processes the file is the reader: it takes
 * records from htslib into batches (counting secondary/supplementary records and applying
 * the -p cutoff), and hands them over a bounded queue to W workers. Each worker converts,
 * parses and collects into a collection of its own; the collections are added up after the
 * workers have joined. Every count is an integer sum (or a min/max), so the result does not
 * depend on which worker saw which batch: the output is the same for any -t.
 * ---------------------------------------------------------------------------------- */

/* Everything one thread needs to turn batches into counts. */
struct file_worker {
        struct batch_queue *q;          /* NULL: the batches are processed by the caller */
        struct tl_seq_buffer *sb;
        struct stat_collection *s;
        int own_s;                      /* the worker allocated s and frees it */
        int is_fastx;
        int ok;
};

/* A sequence buffer with its alignment data. */
static int work_buffer_alloc(struct tl_seq_buffer **buf, int size, int is_fastx)
{
        struct tl_seq_buffer *sb = NULL;

        RUN(alloc_tl_seq_buffer(&sb, size));
        if(!is_fastx){
                RUN(add_aln_data(sb));
        }
        *buf = sb;
        return OK;
ERROR:
        if(sb){
                free_aln_data(sb);
                free_tl_seq_buffer(sb);
        }
        return FAIL;
}

static void work_buffer_free(struct tl_seq_buffer *sb)
{
        if(sb){
                free_aln_data(sb);
                free_tl_seq_buffer(sb);
        }
}

static void file_worker_free(struct file_worker *w)
{
        work_buffer_free(w->sb);
        w->sb = NULL;
        if(w->own_s){
                stat_collection_free(w->s);
        }
        w->s = NULL;
}

/* <s>: the collection to fill, or NULL for a private one. */
static int file_worker_init(struct file_worker *w, struct samstat_param *p, int is_fastx, struct stat_collection *s)
{
        w->q = NULL;
        w->sb = NULL;
        w->s = s;
        w->own_s = 0;
        w->is_fastx = is_fastx;
        w->ok = 1;
        if(!w->s){
                w->own_s = 1;
                RUN(stat_collection_alloc(&w->s));
        }
        RUN(work_buffer_alloc(&w->sb, p->buffer_size, is_fastx));
        return OK;
ERROR:
        file_worker_free(w);
        return FAIL;
}

/* Convert, parse and collect one batch. Touches only the worker's own state. */
static int process_batch(struct file_worker *w, const struct bam_batch *b)
{
        struct tl_seq_buffer *sb = w->sb;

        RUN(convert_bam_batch(b, w->is_fastx, sb));
        sb->L = TL_SEQ_BUFFER_DNA;
        sb->base_quality_offset = 33; /* Can I get this from sam - don't need to part of the standard */
        if(!w->is_fastx){
                /* parse_alignment fails only when it cannot allocate the alignment strings,
                 * and then leaves them NULL: collecting such a record would dereference them,
                 * so the whole batch (and with it the input) fails */
                for(int i = 0; i < sb->num_seq;i++){
                        RUN(parse_alignment(sb->sequences[i]));
                        RUN(fix_orientation(sb->sequences[i]));
                }
        }
        RUN(collect_stats(sb, w->s));
        RUN(trim_seq_buffer(sb, sb->num_seq));
        if(!w->is_fastx){
                clear_aln_data(sb);
        }
        reset_tl_seq_buffer(sb);
        return OK;
ERROR:
        return FAIL;
}

/* A counter on standard error every ~1 M reads (p->progress_step), only for a person
 * watching: standard error must be a terminal. With -v the log lines are written instead. */
static void log_progress(struct samstat_param *p, const char *name, uint64_t n_read, uint64_t *old_n_read)
{
        if(n_read < *old_n_read + p->progress_step){
                return;
        }
        if(p->verbose){
                LOG_MSG("Processed %"PRIu64" sequences", n_read);
        }else if(isatty(STDERR_FILENO)){
                fprintf(stderr, "\rsamstat: %s: %"PRIu64" reads processed", name, n_read);
        }
        *old_n_read = n_read;
}

/* Remove the progress counter again (it is not followed by a newline). */
static void clear_progress(const struct samstat_param *p, uint64_t old_n_read)
{
        if(old_n_read && !p->verbose && isatty(STDERR_FILENO)){
                fputs("\r\033[K", stderr);
        }
}

static double seconds_since(const struct timespec *t0)
{
        struct timespec t1;

        clock_gettime(CLOCK_MONOTONIC, &t1);
        return (double) (t1.tv_sec - t0->tv_sec) + (double) (t1.tv_nsec - t0->tv_nsec) * 1e-9;
}

static const char *format_name(const htsFormat *f)
{
        switch(f->format){
        case sam:
                return "SAM";
        case bam:
                return "BAM";
        case cram:
                return "CRAM";
        case fastq_format:
                return "FASTQ";
        case fasta_format:
                return "FASTA";
        default:
                return NULL;
        }
}

/* "12.4M reads" style counts and "18.2 kb" lengths for the one-line summary. */
static void fmt_count(char *buf, size_t n, uint64_t v)
{
        if(v < 1000){
                snprintf(buf, n, "%"PRIu64, v);
        }else if(v < 1000000){
                snprintf(buf, n, "%.1fk", (double) v / 1e3);
        }else if(v < 1000000000){
                snprintf(buf, n, "%.1fM", (double) v / 1e6);
        }else{
                snprintf(buf, n, "%.1fG", (double) v / 1e9);
        }
}

static void fmt_length(char *buf, size_t n, uint64_t v)
{
        if(v < 1000){
                snprintf(buf, n, "%"PRIu64" bp", v);
        }else if(v < 1000000){
                snprintf(buf, n, "%.1f kb", (double) v / 1e3);
        }else{
                snprintf(buf, n, "%.1f Mb", (double) v / 1e6);
        }
}

/* standard-plots.md section 6: one line per input, written after the report. <s> is
 * finalised by stat_report. One fprintf, so lines of concurrent inputs do not mix. */
static int print_summary_line(const struct samstat_param *p, int id, const struct stat_collection *s, const char *format, uint64_t n_read, double secs)
{
        char *base = NULL;
        char out[2048];
        char reads[32];
        char n50[32];
        char extra[128];
        char partial[96];
        const char *stem = samstat_output_stem(p, id);
        size_t pos = 0;

        if(p->outdir){
                RUN(tld_get_filename((char *) stem, &base));
                snprintf(out, sizeof(out), "%s/%s.samstat.html", p->outdir, base);
        }else{
                snprintf(out, sizeof(out), "%s.samstat.html", stem);
        }
        fmt_count(reads, sizeof(reads), s->n_total);
        extra[0] = 0;
        if(s->n_total){
                if(s->total_bases){
                        fmt_length(n50, sizeof(n50), s->n50_best);
                        pos += (size_t) snprintf(extra + pos, sizeof(extra) - pos, ", N50 %s", n50);
                }
                if(!isnan(s->mean_q)){
                        pos += (size_t) snprintf(extra + pos, sizeof(extra) - pos, ", mean Q %.1f", s->mean_q);
                }
                if(s->is_aligned && !isnan(s->pct_mapped)){
                        snprintf(extra + pos, sizeof(extra) - pos, ", %.1f%% mapped", s->pct_mapped);
                }
        }
        partial[0] = 0;
        if(s->is_partial_report){
                snprintf(partial, sizeof(partial), ", partial: first %"PRIu64" reads", n_read);
        }
        fprintf(stderr, "samstat: %s -> %s (%s%s%s reads%s, %s layout%s, %.0f s)\n",
                stem, out, format[0] ? format : "", format[0] ? ", " : "", reads, extra,
                s->layout == COLLECT_LAYOUT_LONG ? "long" : "short", partial, secs);
        if(base){
                MFREE(base);
        }
        return OK;
ERROR:
        if(base){
                MFREE(base);
        }
        return FAIL;
}

/* The calling thread does everything: -t 1, one input per thread, or no thread could be
 * started. <w> collects directly. */
static int run_inline(struct samstat_param *p, struct sam_bam_file *f, struct file_worker *w, uint64_t *n_read)
{
        struct bam_batch *b = NULL;
        uint64_t old_n_read = 0;

        RUN(bam_batch_alloc(&b, p->buffer_size));
        while(1){
                bam_batch_reset(b);
                RUN(read_bam_batch(f, b));
                if(b->n == 0){
                        break;
                }
                RUN(process_batch(w, b));
                *n_read += (uint64_t) b->n;
                log_progress(p, f->filename, *n_read, &old_n_read);
        }
        bam_batch_free(b);
        return OK;
ERROR:
        bam_batch_free(b);
        return FAIL;
}

#if HAVE_PTHREADS
static void *file_worker_main(void *data)
{
        struct file_worker *w = (struct file_worker *) data;
        struct bam_batch *b = NULL;

        while((b = (struct bam_batch *) batch_queue_get_filled(w->q)) != NULL){
                if(process_batch(w, b) != OK){
                        w->ok = 0;
                        batch_queue_abort(w->q);
                        break;
                }
                bam_batch_reset(b);
                batch_queue_put_free(w->q, b);
        }
        return NULL;
}

/* The calling thread reads, <n_workers> threads process. Their collections are added to
 * <s> after they have joined. */
static int run_pipeline(struct samstat_param *p, struct sam_bam_file *f, struct stat_collection *s, int n_workers, uint64_t *n_read)
{
        struct file_worker *w = NULL;
        struct bam_batch **batches = NULL;
        pthread_t *threads = NULL;
        struct batch_queue *q = NULL;
        uint64_t old_n_read = 0;
        int n_batches = n_workers + 2;  /* one being filled, one waiting, W being processed */
        int n_workers_init = 0;
        int n_batches_init = 0;
        int n_created = 0;
        int reader_ok = 1;
        int ret = FAIL;

        MMALLOC(w, sizeof(struct file_worker) * (size_t) n_workers);
        MMALLOC(threads, sizeof(pthread_t) * (size_t) n_workers);
        MMALLOC(batches, sizeof(struct bam_batch *) * (size_t) n_batches);
        for(n_workers_init = 0; n_workers_init < n_workers; n_workers_init++){
                RUN(file_worker_init(&w[n_workers_init], p, f->is_fastx, NULL));
        }
        for(n_batches_init = 0; n_batches_init < n_batches; n_batches_init++){
                RUN(bam_batch_alloc(&batches[n_batches_init], p->buffer_size));
        }
        RUN(batch_queue_alloc(&q, n_batches, (void **) batches));

        for(int i = 0; i < n_workers; i++){
                /* the workers that run are w[0 .. n_created - 1] */
                int rc;

                w[n_created].q = q;
                rc = pthread_create(&threads[n_created], NULL, file_worker_main, &w[n_created]);
                if(rc){
                        WARNING_MSG("pthread_create() failed with return code %d.", rc);
                        w[n_created].q = NULL;
                }else{
                        n_created++;
                }
        }
        if(n_created == 0){
                /* no thread at all: do the work here */
                RUN(run_inline(p, f, &w[0], n_read));
                RUN(stat_collection_merge(s, w[0].s));
                ret = OK;
                goto ERROR;
        }

        while(1){
                struct bam_batch *b = (struct bam_batch *) batch_queue_get_free(q);

                if(b == NULL){
                        break;          /* a worker failed */
                }
                if(read_bam_batch(f, b) != OK){
                        reader_ok = 0;
                        batch_queue_abort(q);
                        break;
                }
                if(b->n == 0){
                        break;          /* end of the input (or of the -p cutoff) */
                }
                *n_read += (uint64_t) b->n;
                log_progress(p, f->filename, *n_read, &old_n_read);
                batch_queue_put_filled(q, b);
        }
        batch_queue_finish(q);
        for(int i = 0; i < n_created; i++){
                pthread_join(threads[i], NULL);
        }
        if(reader_ok && !batch_queue_aborted(q)){
                for(int i = 0; i < n_created; i++){
                        RUN(stat_collection_merge(s, w[i].s));
                }
                ret = OK;
        }
ERROR:
        /* shared by the success and the error path: every thread has been joined */
        batch_queue_free(q);
        for(int i = 0; i < n_batches_init; i++){
                bam_batch_free(batches[i]);
        }
        for(int i = 0; i < n_workers_init; i++){
                file_worker_free(&w[i]);
        }
        if(batches){
                MFREE(batches);
        }
        if(threads){
                MFREE(threads);
        }
        if(w){
                MFREE(w);
        }
        return ret;
}
#endif

int process_reads_file(struct samstat_param* p, int id)
{
        struct sam_bam_file* f_handle = NULL;
        struct stat_collection* s = NULL;
        struct reader_opts opts;
        struct file_worker inline_worker;
        uint64_t n_read = 0;
        int n_workers = 0;
        int is_stdin = (strcmp(p->infile[id], SAMSTAT_STDIN) == 0);
        struct timespec t0;
        const char *fname;
        char format[16];
        htsFormat hts_fmt;
        int have_hts_fmt = 0;
        int have_inline_worker = 0;
        int ret = FAIL;

        /* readable, not only present: an unreadable file gets the reason here instead of
         * a less specific message from htslib */
        if(!is_stdin && access(p->infile[id], R_OK) != 0){
                ERROR_MSG("cannot read %s: %s", p->infile[id], strerror(errno));
        }

#if HAVE_PTHREADS
        n_workers = samstat_worker_threads(p);
#endif
        opts.pool_threads = samstat_pool_threads(p);
        opts.batch_records = (uint64_t) p->buffer_size;
        opts.batch_bases = p->batch_bases;
        opts.max_primary = p->top_n;
        opts.max_read_len = p->max_read_len;

        clock_gettime(CLOCK_MONOTONIC, &t0);
        RUN(open_reads(&f_handle, p->infile[id], &opts));
        format[0] = 0;
        if(f_handle->in){
                /* copied: the handle is closed before the report is written */
                hts_fmt = *hts_get_format(f_handle->in);
                have_hts_fmt = 1;
        }
        fname = have_hts_fmt ? format_name(&hts_fmt) : NULL;
        if(fname){
                snprintf(format, sizeof(format), "%s", fname);
        }
        if(p->verbose){
                LOG_MSG("%s: workers: %d, decompression threads: %d", p->infile[id], MACRO_MAX(n_workers, 1), opts.pool_threads);
        }

        RUN(stat_collection_alloc(&s));

        /* is_aligned is derived in stat_collection_finalise from has_sq / mapped reads */
        s->has_sq = f_handle->has_sq;

        if(n_workers > 0){
#if HAVE_PTHREADS
                RUN(run_pipeline(p, f_handle, s, n_workers, &n_read));
#endif
        }else{
                /* the one collection is filled directly */
                RUN(file_worker_init(&inline_worker, p, f_handle->is_fastx, s));
                have_inline_worker = 1;
                RUN(run_inline(p, f_handle, &inline_worker, &n_read));
        }
        s->n_secondary = f_handle->n_secondary;
        s->n_supplementary = f_handle->n_supplementary;
        s->is_partial_report = f_handle->partial;
        clear_progress(p, n_read);
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
        RUN(stat_report(s, p, id, have_hts_fmt ? &hts_fmt : NULL));
        RUN(print_summary_line(p, id, s, format, n_read, seconds_since(&t0)));
        ret = OK;
ERROR:
        /* shared by the success and the error path */
        if(f_handle){
                close_bam(f_handle);
        }
        if(have_inline_worker){
                file_worker_free(&inline_worker);
        }
        stat_collection_free(s);
        return ret;
}
