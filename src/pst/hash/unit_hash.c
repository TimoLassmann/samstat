#include "rng/tld-rng.h"
#include "tld.h"
#include "lhash.h"
int main(void)
{

        struct rng_state* rng = NULL;
        void *h = NULL;
        uint8_t* program = NULL;
        int len;
        int i;
        int j;

        int p_len = 10;
        RUN(init_lhash(&h));
        RUN(init_rng(&rng, 0));


        galloc(&program,p_len);


        for (i = 0; i < 10000; ++i) {

                len = tl_random_int(rng, p_len)+1;
                len = MACRO_MAX(0, len);
                len = MACRO_MIN(p_len, len);
                for(j = 0; j < len;j++){
                        program[j] = tl_random_int(rng, 4);
                }

                insert_lhash(h, program, len);
        }

        for (i = 0; i < 10; ++i) {

                len = tl_random_int(rng, p_len)+1;
                len = MACRO_MAX(0, len);
                len = MACRO_MIN(3, len);
                for(j = 0; j < len;j++){
                        program[j] = tl_random_int(rng, 4);
                }
                delete_lhash(h, program, len);
        }

        print_lhash(h);
        free_lhash(h);
        free_rng(rng);
        gfree(program);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
