#include "core/tld-core.h"
#include "seq/tld-seq.h"
#include "string/str.h"
#include "tld.h"

#include "lhash.h"
#include "pst_alloc.h"

#include "pst_structs.h"
#define PST_IMPORT
#include "pst.h"


static int reset_hash(struct pst_model*m);
static int add_to_hash(struct pst_model*m, uint8_t*p, uint32_t l,double w);
static int delete_from_hash(struct pst_model*m, uint8_t*p, uint32_t l, double w);

static int run_build_pst( struct pst_model* model);
static struct pst_node* build_pst(struct pst_model* m,struct pst_node* n);
static struct pst_node* make_flat_pst(struct pst* pst, int maxlen,struct fpst*f,int curf,struct pst_node* n);

static int alloc_pst_model(struct pst_model** model, int L,double min_error, double gamma);
static int reset_pst_model(struct pst_model* m);
static void free_pst_model(struct pst_model* m);

int pst_model_create(struct pst_model **model, struct tl_seq_buffer *sb)
{
        struct pst_model* m = NULL;
        /* struct alphabet* a = NULL; */
        LOG_MSG("L: %d", sb->L);



        RUN(alloc_pst_model(&m, sb->L, -1.0, -1.0));

        for(int i = 0; i < 10;i++){
                for(int j = 0; j < TLD_STRLEN(sb->sequences[i]->seq); j++){
                        fprintf(stdout,"%c",sb->sequences[i]->seq->str[j]);
                }
                fprintf(stdout,"\n");
                convert_to_internal(sb->data, sb->sequences[i]->seq->str, sb->sequences[i]->seq->len);
                for(int j = 0; j < TLD_STRLEN(sb->sequences[i]->seq); j++){
                        fprintf(stdout,"%d",sb->sequences[i]->seq->str[j]);
                }
                fprintf(stdout,"\n");
        }

        exit(0);


        /* RUN(add_to_hash(m, p, l,1)); */

        m->pst->num_nodes = 0;
        RUN(run_build_pst(m));

        *model = m;
        return OK;
ERROR:
        return FAIL;
}

int reset_hash(struct pst_model*m)
{
        uint32_t i;
        for(i = 0;i <= m->depth;i++){
                m->counts_l[i] = 0.0;
        }

        RUN(reset_lhash(m->h));
        return OK;
ERROR:
        return FAIL;
}

int add_to_hash(struct pst_model*m, uint8_t*p, uint32_t l,double w)
{
        uint32_t i;
        uint32_t j;
        ASSERT(m != NULL, "No model");

        if(!m->h){
                RUN(init_lhash(&m->h));
                m->counts_l = NULL;
                galloc(&m->counts_l, m->depth+1);
                for(i = 0;i <= m->depth;i++){
                        m->counts_l[i] = 0.0;
                }
        }
        for(i = 1;i <= m->depth;i++){
                if(l >= i){
                        m->counts_l[i] += l - (i-1);
                }
        }

        for(i = 0; i < l;i++){
                for(j = 1; j < m->depth;j++){
                        if(i +j > l){
                                break;
                        }
                        RUN(insert_lhash(m->h, p+i, j,w));
                }
        }
        return OK;
ERROR:
        return FAIL;
}


int delete_from_hash(struct pst_model*m, uint8_t*p, uint32_t l, double w)
{
        uint32_t i;
        uint32_t j;
        ASSERT(m != NULL, "No model");

        if(!m->h){
                WARNING_MSG("There is no hash to delete items from.");
                return OK;
        }
        for(i = 1;i <= m->depth;i++){
                if(l >= i){
                        m->counts_l[i] -= l - (i-1);
                }
        }

        for(i = 0; i < l;i++){
                for(j = 1; j < m->depth;j++){
                        if(i +j > l){
                                break;
                        }
                        RUN(delete_lhash(m->h, p+i ,j,w));
                }
        }

        return OK;
ERROR:
        return FAIL;
}


int run_build_pst( struct pst_model* model)
{
        struct pst* p = NULL;
        struct pst_node* helper = NULL;
        double sum;
        int i;

        int x;
        uint8_t* suf = NULL;
        uint32_t suf_len;
        uint32_t counts;
        double v;

        galloc(&suf, 1);

        p = model->pst;



        helper = p->root;

        RUN(search_hash(model->h, NULL, 0, &counts,&v));

        sum = 0.0;
        for(i = 0;i < p->L;i++){
                suf[0] = i;
                suf_len = 1;
                RUN(search_hash(model->h, suf, suf_len, &counts,&v));
                /* LOG_MSG("%d %d",i, p->L); */
                helper->probability[i] = (double) counts;
                if(counts){
                        helper->value[i] = v;
                }else{
                        helper->value[i] = 0; /* I have no idea whether this is a winning or losing move */
                }

                sum += counts;
        }
        for(i = 0;i < p->L;i++){
                helper->probability[i] /= sum;
                /* pseudocounts */
                helper->probability[i] = (1.0 - p->gamma_min) * helper->probability[i] + (p->gamma_min * p->background[i]);
                x = p->fpst_root->l;
                p->fpst_root->prob[x][i] = helper->probability[i];
                p->fpst_root->value[x][i] = helper->value[i];
                /* LOG_MSG("%d v:%f p:%f", i, helper->value[i], helper->probability[i]); */
                p->fpst_root->links[x][i] = 0;
        }
        /* exit(0); */

        helper = build_pst(model, helper);
        p->fpst_root->l = 0;
        helper = make_flat_pst(p, model->depth, p->fpst_root, 0, helper);
        p->fpst_root->l++;

        /* print_pst(p, helper); */

        /* exit(0); */
        /* free_pst_node(helper,p->L); */
        /* helper = NULL; */

        /* RUN(prob2scaledprob_fpst(p->fpst_root,p->background,p->L)); */
        /* exit(0); */
        /* *pst = p; */
        model->pst->root = helper;
        gfree(suf);
        return OK;
ERROR:
        gfree(suf);
        return FAIL;
}

struct pst_node* build_pst(struct pst_model* m,struct pst_node* n)
{
        struct pst* pst = NULL;

        uint8_t* suf = NULL;
        uint32_t suf_len;
        uint32_t counts;


        int i;
        int j;
        //int add;

        int len = n->label_len;
        double sum = 0.0;
        double* tmp_counts_s = NULL;
        double* tmp_value_s = NULL;
        double P_ask;
        double P_as;
        double Err;
        double v;

        int maxlen;

        pst = m->pst;

        galloc(&suf, m->depth);
        galloc(&tmp_counts_s, m->L);
        galloc(&tmp_value_s, m->L);

        maxlen = m->depth;
        if(len + 1 < maxlen ){
                for(i = 0; i < pst->L;i++){
                        /* Let's calculate whether extending the suffix by letter A (indexed by i)
                           (i.e. adding it to the left of the existing label)
                           adds statistically significantly information to the PST.

                           p                           What I need is the probability of:
                           1) - P(ASK) probability of observing string ASK;
                           where 'A' is the letter by which I extend the suffix,
                           S is the existing suffix and K is any letter in the alphabet L
                           2) - P(SK) probability of observing string SK
                           3) - P(K | S) - probability of observing K after S

                           Stastistical error = (sum over all K) P(ASK) * log( P(ASK) / ( P(K|S)  * P (AS))
                        */
                        // here I start constructing the ASK by adding A and then copying the existing suffix S
                        suf[0] = i;
                        for(j = 1; j < len+1;j++){
                                suf[j] = n->label[j-1];
                        }
                        suf_len = len+1;

                        RUN(search_hash(m->h, suf, suf_len, &counts,&v));
                        if(counts){
                                P_as = (double) (counts) / (m->counts_l[suf_len]);
                                /* P_as = (float) (c+1) / (float)(h->counts_l[len+1] + pst->L); */
                                Err = 0.0;
                                sum = 0.0;
                                //to calculate the probabilities I add K
                                for(j = 0; j < pst->L;j++){
                                        /* ask_code = ask_code << 5ULL | (uint64_t)j; */
                                        suf[len+1] = j;
                                        suf_len = len+2;
                                        /* now the ask code is complete - I can start looking up their counts to get to P(ask) and P(as) */
                                        RUN(search_hash(m->h, suf, suf_len, &counts,&v));
                                        tmp_counts_s[j] = (double) counts;
                                        if(!counts){
                                                v = 0.0;
                                        }
                                        tmp_value_s[j] = v;
                                        sum+= tmp_counts_s[j];
                                        P_ask = (double) (counts) / (m->counts_l[suf_len]);
                                        if(P_ask != 0.0){

                                                Err += P_ask * log(P_ask /  ( n->probability[j] * P_as));
                                        }


                                }

                                /* Shall we add this suffix ?  */
                                /* NOTE: the check for sum in necessary as it is possible for a suffic P(as)
                                   to occur at the end of sequences. In this case we have P(as) but no string
                                   in the input that is "AS*" */
                                /* LOG_MSG("Sum %f Err %f min%f", sum , Err, pst->p_min); */

                                if (sum > 0.0 && Err > pst->p_min){
                                        for(j = 0; j < pst->L;j++){
                                                tmp_counts_s[j] = tmp_counts_s[j] / sum;
                                                tmp_counts_s[j] = tmp_counts_s[j] * ( 1.0 -  pst->gamma_min) + pst->background[j]* pst->gamma_min;

                                        }
                                        if(!n->next[i]){
                                                RUN(alloc_node(&n->next[i] ,suf,len+1,pst->L));

                                        }
                                        for(j = 0; j < pst->L;j++){
                                                n->next[i]->probability[j] = tmp_counts_s[j];
                                                n->next[i]->value[j] = tmp_value_s[j];
                                        }
                                        pst->num_nodes++;
                                        n->next[i] = build_pst(m, n->next[i]);
                                }
                        }
                }
        }
        gfree(tmp_counts_s);
        gfree(tmp_value_s);
        gfree(suf);
        return n;
ERROR:
        return NULL;
}

struct pst_node* make_flat_pst(struct pst* pst, int maxlen,struct fpst*f,int curf,struct pst_node* n)
{

        int i;
        int j;


        int len;
        len = n->label_len;

        if(len + 1 < maxlen ){
                for(i = 0; i < pst->L;i++){

                        if(n->next[i]){

                                f->links[curf][i] = f->l+1;

                                f->l++;
                                if(f->l== f->m){
                                        resize_fpst(f,pst->L);
                                }

                                for(j = 0; j < pst->L;j++){
                                        f->prob[f->l][j] = n->next[i]->probability[j];
                                        f->value[f->l][j] = n->next[i]->value[j];
                                        f->links[f->l][j] = 0;

                                }
                                n->next[i] = make_flat_pst(pst, maxlen,f, f->links[curf][i], n->next[i]);
                        }

                }
        }
        return n;
ERROR:
        return NULL;
}


int alloc_pst_model(struct pst_model** model, int L,double min_error, double gamma)
{
        struct pst_model* m = NULL;
        MMALLOC(m, sizeof(struct pst_model));
        m->h = NULL;
        m->counts_l = NULL;
        m->L = L;
        m->depth = MAX_PST_MODEL_DEPTH;
        m->pst = NULL;
        m->gamma = 0.1;
        m->min_error = -1.0;    /* this will be set in init pst - now nice! */
        m->gamma = -1.0;

        RUN(init_pst(&m->pst, min_error, gamma, m->depth, L));


        *model = m;
        return OK;
ERROR:
        return FAIL;
}

int reset_pst_model(struct pst_model* m)
{
        RUN(reset_pst(m->pst));
        return OK;
ERROR:
        return FAIL;
}

void free_pst_model(struct pst_model* m)
{
        if(m){
                if(m->h){
                        free_lhash(m->h);
                }
                if(m->pst){
                        free_pst(m->pst);
                }

                gfree(m->counts_l);
                MFREE(m);
        }
}
