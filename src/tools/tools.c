#include "tld.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <inttypes.h>


#define  TOOLS_IMPORT
#include "tools.h"


int sample_wo_replacement(int pop_size, int sample_size,struct rng_state* rng, int* sel)
{
        int n = sample_size;
        int N = pop_size;

        int t = 0; // total input records dealt with
        int m = 0; // number of items selected so far
        int c = 0;
        double r;
        for(int i = 0; i < sample_size;i++){
                sel[i] = 0;
        }
        while (m < n){
                r = tl_random_double(rng);
                if ((N - t)*r >= (n - m) ){
                        t++;
                }else{
                        sel[c] = t;
                        c++;
                        /* samples[m] = t; */
                        t++;
                        m++;
                }
        }
        return OK;
}


int get_file_size(char *filename, uint64_t *size)
{
        uint64_t s = 0;
        ASSERT(filename != NULL, "No filename given");

        ASSERT(tld_file_exists(filename) == OK,"File: %s does not exists",filename);

        struct stat buf;
        int local_ret = 0;
        local_ret = stat(filename, &buf );

        if ( local_ret != 0 ){
                ERROR_MSG("This should never happen..");
        }

        s = (uint64_t)buf.st_size;
        *size =s;

        return OK;
ERROR:
        return FAIL;
}
