/* Unit test for the test-only limits read from the environment (param.c).
 *
 * SAMSTAT_BATCH_RECORDS, SAMSTAT_BATCH_BASES, SAMSTAT_MAX_READ_LEN and SAMSTAT_PROGRESS_STEP
 * may lower the built-in values, never raise them: memory is allocated in proportion to the
 * batch limits, so a large value must be ignored. Only param_init() runs here; nothing is
 * allocated by the values.
 */
#include "tldlite.h"
#include "../../src/param/param.h"
#include <stdio.h>
#include <stdlib.h>

static int n_checks = 0;
static int n_failed = 0;

#define CHECK(COND, ...) do{                                            \
                n_checks++;                                             \
                if(!(COND)){                                            \
                        n_failed++;                                     \
                        fprintf(stderr, "FAIL (%s:%d): %s: ", __FILE__, __LINE__, #COND); \
                        fprintf(stderr, __VA_ARGS__);                   \
                        fprintf(stderr, "\n");                          \
                }                                                       \
        }while(0)

static const char* const vars[4] = {"SAMSTAT_BATCH_RECORDS", "SAMSTAT_BATCH_BASES", "SAMSTAT_MAX_READ_LEN", "SAMSTAT_PROGRESS_STEP"};

static void set_all(const char* v)
{
        for(int i = 0; i < 4; i++){
                if(v){
                        setenv(vars[i], v, 1);
                }else{
                        unsetenv(vars[i]);
                }
        }
}

/* the four values of a fresh parameter set */
static int limits(uint64_t out[4])
{
        struct samstat_param* p = NULL;
        RUN(param_init(&p));
        out[0] = (uint64_t) p->buffer_size;
        out[1] = p->batch_bases;
        out[2] = p->max_read_len;
        out[3] = p->progress_step;
        param_free(p);
        return OK;
ERROR:
        return FAIL;
}

int main(void)
{
        const uint64_t builtin[4] = {SAMSTAT_BATCH_RECORDS, SAMSTAT_BATCH_BASES, SAMSTAT_MAX_READ_LEN, SAMSTAT_PROGRESS_STEP};
        /* values that must leave the built-in limits unchanged */
        const char* const ignored[] = {"99999999999999", "18446744073709551615", "18446744073709551616",
                                       "4294967296", "2147483647", "0", "-5", "abc", "12abc", " 7", ""};
        uint64_t v[4];

        set_all(NULL);
        RUN(limits(v));
        for(int i = 0; i < 4; i++){
                CHECK(v[i] == builtin[i], "unset %s: %llu, expected %llu", vars[i], (unsigned long long) v[i], (unsigned long long) builtin[i]);
        }
        for(size_t k = 0; k < sizeof(ignored) / sizeof(ignored[0]); k++){
                set_all(ignored[k]);
                RUN(limits(v));
                for(int i = 0; i < 4; i++){
                        CHECK(v[i] == builtin[i], "%s='%s': %llu, expected the built-in %llu", vars[i], ignored[k],
                              (unsigned long long) v[i], (unsigned long long) builtin[i]);
                }
        }
        /* lowering works */
        set_all("7");
        RUN(limits(v));
        for(int i = 0; i < 4; i++){
                CHECK(v[i] == 7, "%s=7: %llu", vars[i], (unsigned long long) v[i]);
        }
        set_all(NULL);
        fprintf(stdout, "unit_param: %d checks, %d failed\n", n_checks, n_failed);
        return n_failed ? 1 : 0;
ERROR:
        fprintf(stderr, "unit_param: setup error\n");
        return 2;
}
