#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>

#ifdef THREAD_DATA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct thread_data {
        /* The next input file to take. Each thread claims a file with one atomic
         * fetch-and-add, so there is no lock to take or give back (and no error path
         * that could leave one held). */
        atomic_int next_file;
        uint8_t* ok;            /* per input file: 1 if processed successfully */
        int n_threads;
        struct samstat_param* p;

        pthread_t* threads;
        pthread_attr_t attr;    /* the file threads start with it (samstat_thread_attr_init) */
};
struct samstat_param;

/* Stack size of the threads samstat starts. A file thread decodes CRAM itself when
 * there is no htslib pool, and htslib's rANS codecs need more than 2 MiB of stack
 * (htslib gives its own pool threads at least 3 MiB). The default for a secondary thread
 * is 512 KiB on macOS and 8 MiB on most Linux systems: ask for 8 MiB everywhere, the
 * same as a main thread gets. It is address space, reserved, not memory used. */
#define SAMSTAT_THREAD_STACK ((size_t) 8 * 1024 * 1024)

/* <attr> for a samstat thread: at least SAMSTAT_THREAD_STACK of stack. */
EXTERN int samstat_thread_attr_init(pthread_attr_t *attr);

/* Bounded queue between the reader of one input file and its workers (WP4.1). A fixed
 * set of items (batches) moves around a circle: free -> reader fills -> filled ->
 * worker consumes -> free. One mutex, two condition variables; the reader blocks when
 * every item is in use, the workers when nothing is filled, so memory stays bounded. */
struct batch_queue;
/* The queue does not own the <n> items; all of them start out free. */
EXTERN int batch_queue_alloc(struct batch_queue **q, int n, void **items);
EXTERN void batch_queue_free(struct batch_queue *q);
/* Reader: the next free item; blocks. NULL if the queue was aborted. */
EXTERN void *batch_queue_get_free(struct batch_queue *q);
/* Reader: hand a filled item to the workers. */
EXTERN void batch_queue_put_filled(struct batch_queue *q, void *item);
/* Reader: nothing more follows. Workers drain the filled items, then get NULL. */
EXTERN void batch_queue_finish(struct batch_queue *q);
/* Worker: the next filled item in reading order; blocks. NULL once the reader has
 * finished and every filled item is taken, or at once if the queue was aborted. */
EXTERN void *batch_queue_get_filled(struct batch_queue *q);
/* Worker: give an item back. */
EXTERN void batch_queue_put_free(struct batch_queue *q, void *item);
/* Any thread: stop everything (an error). All blocked threads wake up and get NULL. */
EXTERN void batch_queue_abort(struct batch_queue *q);
EXTERN int batch_queue_aborted(struct batch_queue *q);

EXTERN int thread_data_init(struct thread_data **thread_data, struct samstat_param* param);
EXTERN void thread_data_free(struct thread_data *td);

#undef THREAD_DATA_IMPORT
#undef EXTERN


#endif
