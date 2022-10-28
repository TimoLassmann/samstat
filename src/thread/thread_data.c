#include "tld.h"

#include "../param/param.h"

#define THREAD_DATA_IMPORT
#include "thread_data.h"


int thread_data_init(struct thread_data **thread_data, struct samstat_param* param)
{
        struct thread_data* td = NULL;

        MMALLOC(td, sizeof(struct thread_data));


        td->p = NULL;
        pthread_attr_init(&td->attr);
        pthread_mutex_init(&td->lock, NULL);


        td->n_threads = MACRO_MIN(param->n_infile, param->n_threads);

        galloc(&td->active, param->n_infile);
        for(int i = 0; i < param->n_infile; i++){
                td->active[i] = 0;
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
                if(td->active){
                        gfree(td->active);
                }
                if(td->threads){
                        MFREE(td->threads);
                }
                pthread_mutex_destroy(&td->lock);
                pthread_attr_destroy(&td->attr);
                MFREE(td);
        }
}
