#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#include <inttypes.h>
#include <pthread.h>

#ifdef THREAD_DATA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct thread_data {
        uint8_t* active;
        uint8_t* ok;            /* per input file: 1 if processed successfully */
        pthread_mutex_t lock;
        int n_threads;
        struct samstat_param* p;

        pthread_t* threads;
        pthread_attr_t attr;
};
struct samstat_param;

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
