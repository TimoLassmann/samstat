#include "tldlite.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int header(tld_strbuf *b);
static int footer(tld_strbuf *b);

static int nuc_to_internal(tld_strbuf* b);
static int qual_error_table(tld_strbuf* b);
static int length_edge_tables(tld_strbuf* b);

int main(int argc, char *argv[])
{
        tld_strbuf* b = NULL;
        FILE* f_ptr = NULL;
        if(argc != 2){
                ERROR_MSG("This helper program requires exactly one argument - the name of an output file.");
        }

        RUN(tld_strbuf_alloc(&b,1024));
        RUN(header(b));
        RUN(nuc_to_internal(b));
        RUN(qual_error_table(b));
        RUN(length_edge_tables(b));
        RUN(footer(b));
        f_ptr = fopen(argv[1], "w");
        if(f_ptr == NULL){
                ERROR_MSG("Unable to open file %s for writing.", argv[1]);
        }
        fprintf(f_ptr,"%s", TLD_STR(b));

        fclose(f_ptr);

        tld_strbuf_free(b);

        return EXIT_SUCCESS;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        if(b){
                tld_strbuf_free(b);
        }
        return EXIT_FAILURE;
}

int nuc_to_internal(tld_strbuf* b)
{
        char buffer[1024];
        uint8_t arr[256];
        for(int i = 0; i < 256;i++){
                arr[i] = 4;
        }

        arr['A'] = 0;//	Ala	Alanine
        arr['a'] = 0;//	Ala	Alanine
        arr['C'] = 1;//	Ala	Alanine
        arr['c'] = 1;//	Ala	Alanine
        arr['G'] = 2;//	Ala	Alanine
        arr['g'] = 2;//	Ala	Alanine
        arr['T'] = 3;//	Ala	Alanine
        arr['t'] = 3;//	Ala	Alanine


        RUN(tld_append(b, "const uint8_t nuc_to_internal[256] = {\n"));

        for(int i = 0; i < 256;i++){
                snprintf(buffer, 1024,"%3d,", arr[i]);
                RUN(tld_append(b,buffer));
                if(i != 0 && (i +1) % 16 == 0){
                        RUN(tld_append_char(b, '\n'));
                }else  if(i != 0 && (i+1) %4 == 0){
                        RUN(tld_append_char(b, ' '));
                }
        }
        b->len--;
        RUN(tld_append(b, "\n};\n"));
        return OK;
ERROR:
        return FAIL;
}


/* standard-plots.md 3.3: E[q] = round(10^(-q/10) * 2^30) for q = 0..93, so that no pow()
 * runs at collect time. The values are not close to a rounding tie (the closest is 0.015
 * away), so libm differences in the last bit cannot change them. */
#define N_QUAL 94
int qual_error_table(tld_strbuf* b)
{
        char buffer[64];

        RUN(tld_append(b, "#define CT_N_QUAL 94\n"));
        RUN(tld_append(b, "const uint32_t qual_error[CT_N_QUAL] = {\n"));
        for(int q = 0; q < N_QUAL; q++){
                double e = pow(10.0, -(double) q / 10.0) * 1073741824.0;
                snprintf(buffer, 64, "%10lluu,", (unsigned long long) llround(e));
                RUN(tld_append(b, buffer));
                if((q + 1) % 8 == 0){
                        RUN(tld_append_char(b, '\n'));
                }else{
                        RUN(tld_append_char(b, ' '));
                }
        }
        b->len--;
        RUN(tld_append(b, "\n};\n"));
        return OK;
ERROR:
        return FAIL;
}

/* The master grid of 3.1: g(e0) = e0, then e[j+1] = max(e[j] + 1, (e[j] * 101 + 99) / 100)
 * up to and including the first edge >= 10,000,000. Integer arithmetic only. */
#define GRID_STOP 10000000ull

static uint64_t grid_next(uint64_t e)
{
        uint64_t n = (e * 101ull + 99ull) / 100ull;
        return (n > e + 1) ? n : e + 1;
}

/* Storage edges: lengths 0..999, then g(1000). Display edges (plot 1 in the long layout,
 * MultiQC): 0, g(1) below 1000, then the storage edges from 1000 upwards, so that every
 * display bin is the union of whole storage bins. */
int length_edge_tables(tld_strbuf* b)
{
        char buffer[64];
        uint64_t e;
        int n_storage = 1000;
        int n_display = 1;
        int col = 0;

        for(e = 1000; ; e = grid_next(e)){
                n_storage++;
                if(e >= GRID_STOP){
                        break;
                }
        }
        for(e = 1; e < 1000; e = grid_next(e)){
                n_display++;
        }
        n_display += n_storage - 1000;

        snprintf(buffer, 64, "#define CT_N_LENGTH_EDGES %d\n", n_storage);
        RUN(tld_append(b, buffer));
        RUN(tld_append(b, "const uint32_t length_edges[CT_N_LENGTH_EDGES] = {\n"));
        for(int i = 0; i < 1000; i++){
                snprintf(buffer, 64, "%u,", (unsigned) i);
                RUN(tld_append(b, buffer));
                RUN(tld_append_char(b, (++col % 16 == 0) ? '\n' : ' '));
        }
        for(e = 1000; ; e = grid_next(e)){
                snprintf(buffer, 64, "%u,", (unsigned) e);
                RUN(tld_append(b, buffer));
                RUN(tld_append_char(b, (++col % 16 == 0) ? '\n' : ' '));
                if(e >= GRID_STOP){
                        break;
                }
        }
        b->len--;
        RUN(tld_append(b, "\n};\n"));

        col = 0;
        snprintf(buffer, 64, "#define CT_N_DISPLAY_EDGES %d\n", n_display);
        RUN(tld_append(b, buffer));
        RUN(tld_append(b, "const uint32_t display_edges[CT_N_DISPLAY_EDGES] = {\n0,"));
        RUN(tld_append_char(b, ' '));
        for(e = 1; e < 1000; e = grid_next(e)){
                snprintf(buffer, 64, "%u,", (unsigned) e);
                RUN(tld_append(b, buffer));
                RUN(tld_append_char(b, (++col % 16 == 0) ? '\n' : ' '));
        }
        for(e = 1000; ; e = grid_next(e)){
                snprintf(buffer, 64, "%u,", (unsigned) e);
                RUN(tld_append(b, buffer));
                RUN(tld_append_char(b, (++col % 16 == 0) ? '\n' : ' '));
                if(e >= GRID_STOP){
                        break;
                }
        }
        b->len--;
        RUN(tld_append(b, "\n};\n"));
        return OK;
ERROR:
        return FAIL;
}

int header(tld_strbuf *b)

{
        RUN(tld_append(b, "#ifndef CONVERT_TABLES\n"));
        RUN(tld_append(b, "#define CONVERT_TABLES\n"));
        RUN(tld_append(b, "#include <stdint.h>\n"));
        return OK;
ERROR:
        return FAIL;
}
int footer(tld_strbuf *b)
{
        RUN(tld_append(b, "#endif\n"));
        return OK;
ERROR:
        return FAIL;
}
