#include "tld.h"
#include <stdio.h>
#include <stdlib.h>

static int header(tld_strbuf *b);
static int footer(tld_strbuf *b);

static int nuc_to_internal(tld_strbuf* b);

int main(int argc, char *argv[])
{
        tld_strbuf* b = NULL;
        FILE* f_ptr = NULL;
        if(argc != 2){
                ERROR_MSG("This helper program requires exactly one argument - the name of an output file.");
        }

        tld_strbuf_alloc(&b,1024);
        RUN(header(b));
        RUN(nuc_to_internal(b));
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
                        tld_append_char(b, '\n');
                }else  if(i != 0 && (i+1) %4 == 0){
                        tld_append_char(b, ' ');
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
