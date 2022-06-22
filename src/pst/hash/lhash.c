#include "core/tld-core.h"
#include "tld.h"
#include <stdint.h>

#define  LHASH_IMPORT
#include "lhash.h"


#include "khash.h"


#include "murmur3.h"

struct lhash_item{
        uint32_t* fragment;
        uint32_t counts;
        double wins;
        uint32_t len;
};

/* static inline uint32_t get_bitshift_hashkey(const struct lhash_item *a); */
static inline uint32_t get_murmur(const struct lhash_item* a);
static inline int lhash_items_eq (const struct lhash_item* a,const struct lhash_item* b);


static int create_lhash_item(struct lhash_item** item, uint8_t* p, uint32_t len, double w);
static void free_lhash_item(struct lhash_item* a);


typedef struct lhash_item* tl_item_t;
#define KHASH_MAP_INIT_LHASH(name)                             \
        KHASH_INIT(name, tl_item_t, char, 0,  get_murmur, lhash_items_eq)

KHASH_MAP_INIT_LHASH(phash)

int search_hash(void *hash, uint8_t* p, uint32_t len, uint32_t* counts, double*v)
{
        khash_t(phash) *h = (khash_t(phash)*) hash;
        struct lhash_item* item = NULL;
        khint_t k;

        double value;
        double c;
        /* int ret; */

        *counts = 0;
        RUN(create_lhash_item(&item, p, len,0));


        k = kh_get(phash, h, item);
        if( k != kh_end(h)){
                *counts = kh_key(h,k)->counts;
                value = (kh_key(h,k)->wins);
                c = (double) *counts;
                /* *v = ((value  )  - (c - value  )) / (c); */
                *v = value / c;
                /* LOG_MSG("Value: %f %f %f", *v, value ,c); */
        }



        free_lhash_item(item);
        return OK;
ERROR:
        return FAIL;
}

int insert_lhash(void *hash,uint8_t* p, uint32_t len,double w )
{
        khash_t(phash) *h = (khash_t(phash)*) hash;
        struct lhash_item* item = NULL;
        uint32_t* f = NULL;
        khint_t k;

        int ret;

        
        RUN(create_lhash_item(&item, p, len,w));

        /* uint32_t i; */
        /* fprintf(stdout,"INSERTING : "); */
        /* for(i = 0; i < item->len;i++){ */
        /*         fprintf(stdout,"%d ", item->fragment[i]); */
        /* } */
        /* fprintf(stdout,"\n"); */


        k = kh_put(phash, h, item, &ret);
        if(ret){
                kh_key(h, k) = item;
        }else{
                kh_key(h, k)->counts++;
                kh_key(h, k)->wins += item->wins;

                free_lhash_item(item);
        }

        /* print_lhash(h); */
        return OK;
ERROR:
        if(f){
                MFREE(f);
        }
        return FAIL;
}

int delete_lhash(void* hash,uint8_t* p, uint32_t len, double w)
{
        khash_t(phash) *h = (khash_t(phash)*) hash;
        struct lhash_item* item = NULL;

        khiter_t k;
        int ret;

        RUN(create_lhash_item(&item, p, len,w));

        /* uint32_t i; */
        /* fprintf(stdout,"DELETING : "); */
        /* for(i = 0; i < item->len;i++){ */
        /*         fprintf(stdout,"%d ", item->fragment[i]); */
        /* } */
        /* fprintf(stdout,"\n"); */

        k = kh_put(phash, h, item, &ret);
        if(!ret){
                kh_key(h, k)->counts--;
                kh_key(h, k)->wins -= item->wins;

                if(kh_key(h, k)->counts == 0){
                        free_lhash_item(kh_key(h, k));
                        kh_del(phash, h, k);
                }

        }
        free_lhash_item(item);

        return OK;
ERROR:
        free_lhash_item(item);
        return FAIL;
}

int init_lhash(void** hash)
{
        khash_t(phash) *h = NULL;

        h = kh_init(phash);
        if(h == NULL){
                ERROR_MSG("kn_init failed");
        }

        *hash = h;
        return OK;
ERROR:
        return FAIL;
}

void print_lhash(void* hash)
{

        khash_t(phash) *h = (khash_t(phash)*) hash;
        struct lhash_item* item = NULL;
        khiter_t k;
        uint32_t j;

        for (k = kh_begin(h); k != kh_end(h); ++k){
                if (kh_exist(h, k)){

                        item = kh_key(h, k);

                        /* kh_val(h,k).counts */
                        /* fprintf(stdout, "%d %d\n",k,kh_value(h, k)); */
                        /* hi = kh_value(h, k); */
                        fprintf(stdout, "%d counts: %d val: %f : ",k,item->counts,item->wins);
                        for(j = 0; j < item->len;j++){
                                fprintf(stdout,"%d ", item->fragment[j]);
                        }
                        fprintf(stdout,"\n");


                        //kh_value(h, k) = 1;
                }
        }
}

int reset_lhash(void* hash)
{
        khash_t(phash) *h = (khash_t(phash)*) hash;
        if(h){
                khiter_t k;
                for (k = kh_begin(h); k != kh_end(h); ++k){
                        if (kh_exist(h, k)){
                                /* kh_key(h, k)->counts = 0 */
                                free_lhash_item(kh_key(h, k));
                        }
                }
                kh_clear(phash, h);
        }
        return OK;
}

void free_lhash(void* hash)
{

        khash_t(phash) *h = (khash_t(phash)*) hash;
        if(h){
                khiter_t k;
                for (k = kh_begin(h); k != kh_end(h); ++k){
                        if (kh_exist(h, k)){
                                free_lhash_item(kh_key(h, k));
                        }
                }
                kh_destroy(phash, h);
        }
}

int create_lhash_item(struct lhash_item** item, uint8_t* p, uint32_t len, double w)
{
        struct lhash_item* a = NULL;
        uint32_t i;
        MMALLOC(a, sizeof(struct lhash_item));
        a->fragment = NULL;
        /* LOG_MSG("%d len",len); */
        MMALLOC(a->fragment, sizeof(uint32_t) * len);

        a->len = len;
        for(i = 0; i < len;i++){
                /* LOG_MSG("%d %d", i , p[i]); */
                a->fragment[i] = p[i];
        }
        a->counts = 1;
        a->wins = w;
        *item = a;
        return OK;
ERROR:
        free_lhash_item(a);
        return FAIL;
}


void free_lhash_item(struct lhash_item* a)
{
        if(a){
                if(a->fragment){
                        MFREE(a->fragment);
                }
                MFREE(a);
        }
}

/* static inline uint32_t get_bitshift_hashkey(const struct lhash_item *a) */
/* { */
/*         /\* uint32_t key = 0; *\/ */
/*         /\* uint8_t* x = a->fragment; *\/ */
/*         /\* uint32_t len = a->len; *\/ */
/*         /\* for(int i = 0; i < len;i++){ *\/ */
/*         /\*         key = (key << 2) & x[i]; *\/ */
/*         /\* } *\/ */
/*         /\* return key; *\/ */
/* } */

static inline uint32_t get_murmur(const struct lhash_item* a)
{
        uint32_t* x = a->fragment;
        uint32_t len = a->len;
        uint32_t key;
        len = len << 2;


        MurmurHash3_x86_32((void*)x, len, 42, &key);
        return key;

}
static inline int lhash_items_eq (const struct lhash_item* a,const struct lhash_item* b)
{

        if (a->len == b->len){
                int32_t len = a->len;
                len = len << 2;
                return memcmp (a->fragment,b->fragment, len) == 0;
        }

        return 0;
}

#ifdef LHASH_TEST

#include "tlrng.h"
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
        RUNP(rng = init_rng(0));

        galloc(&program,p_len);


        for (i = 0; i < 10000; ++i) {

                len = tl_random_int(rng, p_len)+1;
                len = MACRO_MAX(0, len);
                len = MACRO_MIN(p_len, len);
                for(j = 0; j < len;j++){
                        program[j] = tl_random_int(rng, 4);
                }

                insert_lhash(h, program, len,1);
        }

        for (i = 0; i < 10; ++i) {

                len = tl_random_int(rng, p_len)+1;
                len = MACRO_MAX(0, len);
                len = MACRO_MIN(3, len);
                for(j = 0; j < len;j++){
                        program[j] = tl_random_int(rng, 4);
                }
                delete_lhash(h, program, len,1);
        }

        print_lhash(h);
        free_lhash(h);
        free_rng(rng);
        gfree(program);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

#endif
