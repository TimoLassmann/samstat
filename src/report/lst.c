#include "alloc/tld-alloc.h"
#include "core/tld-core.h"
#include "string/str.h"
#define LST_IMPORT
#include "lst.h"

int lst_append(lst_node **head, char *content)
{

        lst_node* n = *head;
        if(!n){
                RUN(lst_node_alloc(&n));
                RUN(tld_append(n->data, content));
                *head = n;
        }else{
                while(n->next){
                        n = n->next;
                }

                RUN(lst_node_alloc(&n->next));
                n = n->next;
                RUN(tld_append(n->data, content));
        }
        return OK;
ERROR:
        return FAIL;
}


int lst_concatenate(lst_node *n, char* sep)
{
        if(!n){
                ERROR_MSG("Nothing to concatenate");
        }
        lst_node* head = NULL;

        head = n;

        if(sep){
                tld_append(head->data, sep);
        }

        while(n->next){
                n = n->next;
                tld_append(head->data, TLD_STR(n->data));
                if(sep && n->next){
                        tld_append(head->data, sep);
                }
        }
        lst_free(head->next);
        head->next = NULL;
        return OK;
ERROR:
        return FAIL;
}


int lst_print(lst_node *n)
{
        while(n){
                LOG_MSG("%s",n->data->str);
                n = n->next;
        }
        return OK;
}

void lst_free(lst_node *n)
{
        lst_node* tmp = NULL;

        while(n){
                tmp =n;
                n = n->next;
                lst_node_free(tmp);
        }
}

int lst_node_alloc(lst_node **list_n)
{
        lst_node* n  = NULL;
        MMALLOC(n, sizeof(lst_node));

        RUN(tld_strbuf_alloc(&n->data, 128));
        n->next = NULL;

        *list_n = n;
        return OK;
ERROR:
        return FAIL;
}


void lst_node_free(lst_node *n)
{
        if(n){
                if(n->data){
                        tld_strbuf_free(n->data);
                }
                MFREE(n);
        }
}
