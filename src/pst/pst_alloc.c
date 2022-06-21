#include "tld.h"
#include "pst_structs.h"

#define PST_ALLOC_IMPORT
#include "pst_alloc.h"

int init_pst(struct pst** pst, double min_error, double expected_error, int len, int L)
{
        struct pst* p = NULL;
        int x;
        int i;
        MMALLOC(p, sizeof(struct pst));
        p->root = NULL;
        p->max_observed_len = -1;
        p->len = len;
        p->L = L;
        p->gamma_min = 0.02;
        p->num_nodes = 0;
        if(expected_error != -1.0){
                p->gamma_min = expected_error;
        }
        if(p->gamma_min > 1.0){
                WARNING_MSG("expected error is too high.");
                p->gamma_min = 1.0 ;
        }


        p->p_min = 0.0001;

        if(min_error!= -1){
                p->p_min = min_error;
        }

        //p->r = 0.0001f;
        p->background = NULL;

        p->fpst_root = NULL;

        RUN(alloc_fpst(&p->fpst_root, 64,p->L));

        /* RUN(get_null_model_emissions(&tmp, p->L)); */

        galloc(&p->background,p->L);
        for(i = 0; i < p->L;i++){
                p->background[i] = 1.0/ (double) p->L;
        }

        /* Init first node */
        RUN(alloc_node(&p->root,NULL,0,p->L));
        for(i = 0;i < p->L;i++){

                p->root->probability[i] = p->background[i];
                x = p->fpst_root->l;
                p->fpst_root->prob[x][i] = p->background[i];
                p->fpst_root->value[x][i] = 0.0;
                p->fpst_root->links[x][i] = 0;
        }
        /* p->fpst_root->l = 0; */
        /* p->root = make_flat_pst(p, len ,p->fpst_root, 0, p->root); */
        /* p->fpst_root->l++; */
        *pst = p;


        return OK;
ERROR:
        free(p);
        return FAIL;
}

int reset_pst(struct pst* p)
{
        struct pst_node* n =  NULL;
        int i;
        for(i = 0; i < p->L;i++){
                p->background[i] = 1.0/ (double) p->L;
        }
        p->num_nodes = 0;
        p->fpst_root->l = 0;
        n = p->root;
        for(i = 0;i < p->L;i++){
                /* reset f_pst */
                p->fpst_root->prob[0][i] = p->background[i];
                p->fpst_root->value[0][i] = 0.0;
                p->fpst_root->links[0][i] = 0;

                /* reset pst root  */
                if(n->next[i]){
                        free_pst_node(n->next[i], p->L);

                }
                n->next[i] = NULL;
                n->value[i] = 0.0;
                /* Should I really set this here?  */
                n->probability[i] = p->background[i];
        }
        return OK;
}

void free_pst(struct pst* p)
{
        if(p){
                if(p->fpst_root){
                        free_fpst(p->fpst_root);
                }
                if(p->root){
                        free_pst_node(p->root, p->L);
                }
                gfree(p->background);
                MFREE(p);
        }
}

void free_pst_node(struct pst_node* n, int L)
{
        int i;
        if(n){
                for(i = 0;i < L;i++){
                        if(n->next[i]){
                                free_pst_node(n->next[i],L);
                        }
                }
                gfree(n->value);
                gfree(n->probability);
                gfree(n->label);
                MFREE(n->next);
                MFREE(n);
        }
}


int alloc_node(struct pst_node** node,uint8_t* string,int len, int L)
{
        struct pst_node* n =  NULL;
        int i;
        MMALLOC(n, sizeof(struct pst_node));
        n->label = NULL;

        if(!len){
                galloc(&n->label, len+1);
        }else{
                galloc(&n->label, len);
        }
        n->probability = NULL;

        n->value = NULL;
        n->next = NULL;

        for(i = 0; i < len;i++){
                n->label[i] = string[i];
        }
        n->label_len = len;

        galloc(&n->probability, L);
        galloc(&n->value, L);


        MMALLOC(n->next, sizeof(struct pst_node*) * L);

        for(i =0; i < L;i++){
                n->next[i] = NULL;
                n->value[i] = 0.0;
                /* Should I really set this here?  */
                n->probability[i] = 1.0 / (double) L;
        }
        *node = n;
        return OK;
ERROR:
        return FAIL;
}


int alloc_fpst(struct fpst** fast_pst, int size, int L)
{
        struct fpst*f = NULL;

        ASSERT(L > 3, "Strange alphabet...");
        ASSERT(size > 0, "no space!");
        MMALLOC(f, sizeof(struct fpst));
        f->prob = NULL;
        f->value = NULL;
        f->links = NULL;
        f->m = size;
        f->l= 0;

        galloc(&f->prob,f->m,L);
        galloc(&f->value ,f->m,L);
        galloc(&f->links,f->m,L);

        *fast_pst = f;
        return OK;
ERROR:
        return FAIL;
}

int resize_fpst(struct fpst* f, int L)
{
        int i,j, old;

        old = f->m;
        f->m = f->m + f->m /2;
        galloc(&f->prob,f->m,L);
        galloc(&f->value,f->m,L);
        galloc(&f->links,f->m,L);

        for(i = old;i < f->m;i++){
                for(j = 0;j < L;j++){
                        f->prob[i][j] = 0.0;
                        f->value[i][j] = 0.0;
                        f->links[i][j] = 0;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

void free_fpst(struct fpst* f)
{
        if(f){
                gfree(f->prob);
                gfree(f->value);
                gfree(f->links);
                MFREE(f);
        }
}
