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
        pthread_mutex_t lock;
        int n_threads;
        struct samstat_param* p;

        pthread_t* threads;
        pthread_attr_t attr;
};
struct samstat_param;
EXTERN int thread_data_init(struct thread_data **thread_data, struct samstat_param* param);
EXTERN void thread_data_free(struct thread_data *td);

#undef THREAD_DATA_IMPORT
#undef EXTERN


#endif
