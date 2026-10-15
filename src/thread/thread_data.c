#include "tldlite.h"

#include "../param/param.h"

#define THREAD_DATA_IMPORT
#include "thread_data.h"


int samstat_thread_attr_init(pthread_attr_t *attr)
{
        size_t stack_size = 0;

        ASSERT(pthread_attr_init(attr) == 0, "Cannot initialise the thread attributes");
        if(pthread_attr_getstacksize(attr, &stack_size) != 0 || stack_size < SAMSTAT_THREAD_STACK){
                /* a larger default (ulimit -s on Linux) is kept */
                if(pthread_attr_setstacksize(attr, SAMSTAT_THREAD_STACK) != 0){
                        pthread_attr_destroy(attr);
                        ERROR_MSG("Cannot set the thread stack size to %zu bytes", SAMSTAT_THREAD_STACK);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int thread_data_init(struct thread_data **thread_data, struct samstat_param* param)
{
        struct thread_data* td = NULL;

        MMALLOC(td, sizeof(struct thread_data));


        td->p = NULL;
        td->ok = NULL;
        td->threads = NULL;
        atomic_init(&td->next_file, 0);
        if(samstat_thread_attr_init(&td->attr) != OK){
                MFREE(td);
                ERROR_MSG("Cannot prepare the file threads");
        }


        /* T inputs at once (one per thread), or one after another with all T threads on each */
        td->n_threads = samstat_file_threads(param);

        galloc(&td->ok, param->n_infile);
        for(int i = 0; i < param->n_infile; i++){
                td->ok[i] = 0;
        }
        MMALLOC(td->threads, sizeof(pthread_t) * td->n_threads);

        td->p = param;
        *thread_data = td;
        return OK;
ERROR:
        if(td){
                thread_data_free(td);
        }
        return FAIL;
}

void thread_data_free(struct thread_data *td)
{
        if(td){
                if(td->ok){
                        gfree(td->ok);
                }
                if(td->threads){
                        MFREE(td->threads);
                }
                pthread_attr_destroy(&td->attr);
                MFREE(td);
        }
}


struct batch_queue {
        pthread_mutex_t lock;
        pthread_cond_t has_free;      /* the reader waits for a free item or an abort */
        pthread_cond_t has_filled;    /* workers wait for a filled item, the end or an abort */
        void **free_items;            /* stack of free items */
        void **filled_items;          /* ring buffer of filled items, in reading order */
        int n;
        int n_free;
        int head;
        int n_filled;
        int finished;
        int aborted;
};

int batch_queue_alloc(struct batch_queue **q, int n, void **items)
{
        struct batch_queue *b = NULL;

        ASSERT(n > 0, "A queue needs at least one item");
        MMALLOC(b, sizeof(struct batch_queue));
        b->free_items = NULL;
        b->filled_items = NULL;
        b->n = n;
        b->n_free = 0;
        b->head = 0;
        b->n_filled = 0;
        b->finished = 0;
        b->aborted = 0;
        MMALLOC(b->free_items, sizeof(void *) * (size_t) n);
        MMALLOC(b->filled_items, sizeof(void *) * (size_t) n);
        for(int i = 0; i < n; i++){
                b->free_items[b->n_free++] = items[i];
        }
        pthread_mutex_init(&b->lock, NULL);
        pthread_cond_init(&b->has_free, NULL);
        pthread_cond_init(&b->has_filled, NULL);
        *q = b;
        return OK;
ERROR:
        if(b){
                if(b->free_items){
                        MFREE(b->free_items);
                }
                if(b->filled_items){
                        MFREE(b->filled_items);
                }
                MFREE(b);
        }
        return FAIL;
}

void batch_queue_free(struct batch_queue *q)
{
        if(q){
                pthread_cond_destroy(&q->has_free);
                pthread_cond_destroy(&q->has_filled);
                pthread_mutex_destroy(&q->lock);
                MFREE(q->free_items);
                MFREE(q->filled_items);
                MFREE(q);
        }
}

void *batch_queue_get_free(struct batch_queue *q)
{
        void *item = NULL;

        pthread_mutex_lock(&q->lock);
        while(q->n_free == 0 && !q->aborted){
                pthread_cond_wait(&q->has_free, &q->lock);
        }
        if(!q->aborted){
                item = q->free_items[--q->n_free];
        }
        pthread_mutex_unlock(&q->lock);
        return item;
}

void batch_queue_put_filled(struct batch_queue *q, void *item)
{
        pthread_mutex_lock(&q->lock);
        q->filled_items[(q->head + q->n_filled) % q->n] = item;
        q->n_filled++;
        pthread_cond_signal(&q->has_filled);
        pthread_mutex_unlock(&q->lock);
}

void batch_queue_finish(struct batch_queue *q)
{
        pthread_mutex_lock(&q->lock);
        q->finished = 1;
        pthread_cond_broadcast(&q->has_filled);
        pthread_mutex_unlock(&q->lock);
}

void *batch_queue_get_filled(struct batch_queue *q)
{
        void *item = NULL;

        pthread_mutex_lock(&q->lock);
        while(q->n_filled == 0 && !q->finished && !q->aborted){
                pthread_cond_wait(&q->has_filled, &q->lock);
        }
        if(!q->aborted && q->n_filled > 0){
                item = q->filled_items[q->head];
                q->head = (q->head + 1) % q->n;
                q->n_filled--;
        }
        pthread_mutex_unlock(&q->lock);
        return item;
}

void batch_queue_put_free(struct batch_queue *q, void *item)
{
        pthread_mutex_lock(&q->lock);
        q->free_items[q->n_free++] = item;
        pthread_cond_signal(&q->has_free);
        pthread_mutex_unlock(&q->lock);
}

void batch_queue_abort(struct batch_queue *q)
{
        pthread_mutex_lock(&q->lock);
        q->aborted = 1;
        pthread_cond_broadcast(&q->has_free);
        pthread_cond_broadcast(&q->has_filled);
        pthread_mutex_unlock(&q->lock);
}

int batch_queue_aborted(struct batch_queue *q)
{
        int a;

        pthread_mutex_lock(&q->lock);
        a = q->aborted;
        pthread_mutex_unlock(&q->lock);
        return a;
}
