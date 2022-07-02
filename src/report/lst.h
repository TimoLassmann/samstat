#ifndef LST_H
#define LST_H

#include "string/str.h"
#include "tld.h"

#ifdef LST_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

typedef struct lst_node {
        tld_strbuf* data;
        struct lst_node* next;
} lst_node;


EXTERN int lst_concatenate(lst_node *n, char* sep);
EXTERN int lst_append(lst_node **head, char *content);
EXTERN int lst_print(lst_node *n);
EXTERN void lst_free(lst_node *n);

EXTERN int lst_node_alloc(lst_node **list_n);
EXTERN void lst_node_free(lst_node *n);

#undef LST_IMPORT
#undef EXTERN


#endif
