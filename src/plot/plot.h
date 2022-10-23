#ifndef PLOT_H
#define PLOT_H


#include <stdint.h>

typedef struct tld_string_buffer  tld_strbuf;

#ifdef PLOT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct plot_data {
        uint64_t** data;
        uint64_t** clu_data;
        char** series_label;
        char** group_label;
        char** x_axis_labels;
        uint8_t** visible_tf;
        int32_t n_group;
        int len;
        int clu_len;
        int L;
        tld_strbuf* title;
        tld_strbuf* xlabel;
        tld_strbuf* ylabel;
        tld_strbuf* id;
        tld_strbuf* save_file_name;
        uint8_t* is_plot;
        uint8_t* is_visible;
        int32_t bin_start;
        int32_t bin_size;
        int8_t group_size;
        int8_t viz;
        int8_t type;
        int8_t mod;
        uint8_t x_is_categorical;
        uint32_t target_n_clu;
};

#define PLOT_TYPE_SCATTER 0
#define PLOT_TYPE_LINES 1
#define PLOT_TYPE_BAR 2

#define PLOT_MOD_NORMAL 0
#define PLOT_MOD_ERROR_BAR 1
#define PLOT_MOD_DENSITY 2

#define PLOT_VIZ_ALL 0
#define PLOT_VIZ_FIRSTGROUP 1

#define PLOT_TYPE_SHIFT 8
#define PLOT_TYPE_MASK 0xFF

EXTERN int plot_add(tld_strbuf *o, struct plot_data *d);

EXTERN int plot_data_config(struct plot_data *d, int8_t type, int8_t mod, int8_t group_size, int viz , char* id, char* title, char* x_label, char* y_label, char* savename,char** series_lab, char** group_lab );
EXTERN int plot_data_alloc(struct plot_data **plot_data, int x, int y);
EXTERN int plot_data_resize_len(struct plot_data* pd,int x);
EXTERN void plot_data_free(struct plot_data *pd);

EXTERN int html_header(tld_strbuf *out_buffer, char *filename);
EXTERN int html_end(tld_strbuf *out_buffer);

#undef PLOT_IMPORT
#undef EXTERN


#endif
