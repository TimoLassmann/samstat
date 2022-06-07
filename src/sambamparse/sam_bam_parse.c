
#include "alloc/tld-alloc.h"
#include "core/tld-core.h"
#include "string/str.h"
#include "tld.h"
#include "../htsinterface/htsglue.h"
#include <ctype.h>
#include <string.h>

#define SAM_BAM_PARSE_IMPORT
#include "sam_bam_parse.h"


int get_readable_cigar(struct aln_data *a, char **cigar)
{
        char buffer[1024];
        char* ret = NULL;
        tld_strbuf* b = NULL;
        tld_strbuf_alloc(&b,16);

        char Op;
        int Oplen;
        for(int i = 0;i < a->n_cigar;i++){
                Op = bam_cigar_opchr(a->cigar[i]);
                Oplen = bam_cigar_oplen(a->cigar[i]);
                snprintf(buffer,1024,"%d%c", Oplen,Op);
                tld_append(b, buffer);
        }


        /* LOG_MSG("CIGAR: %s", TLD_STR(b)); */

        galloc(&ret,b->len+1);

        memcpy(ret, TLD_STR(b),TLD_STRLEN(b));
        ret[b->len] = 0;
        /* LOG_MSG("CIGAR: %s", ret); */

        *cigar = ret;

        tld_strbuf_free(b);

        return OK;
ERROR:
        tld_strbuf_free(b);
        gfree(ret);
        return FAIL;
}

int parse_alignment(struct tl_seq *s)
{
        struct aln_data* a = NULL;
        uint8_t* genome = NULL;
        uint8_t* read = NULL;
        int sp;
        int rp;
        char Op;
        int Oplen;
        int aln_len;

        a = s->data;

        aln_len = 0;

        for(int i = 0;i < a->n_cigar;i++){
                Oplen = bam_cigar_oplen(a->cigar[i]);
                aln_len += Oplen;
        }
        aln_len++;


        galloc(&genome, aln_len);
        galloc(&read, aln_len);
        for(int i = 0;i < aln_len;i++){
                genome[i] = 0;
                read[i] = 0;
        }

        rp = 0;
        sp = 0;
        for(int i = 0;i < a->n_cigar;i++){
                Op = bam_cigar_opchr(a->cigar[i]);
                Oplen = bam_cigar_oplen(a->cigar[i]);
                switch (Op) {
                case 'N':
                        break;
                case 'P':
                        /* We'll ignore padding */
                        break;
                case 'H':
                        /* We'll ignore hard clipping as we have nothing to work on */
                        break;
                case 'S':
                        sp++;
                        break;
                case 'M':
                case '=':
                case 'X':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = s->seq[sp];
                                sp++;
                                rp++;
                        }
                        /* tld_append_char(genome, s->seq[sp]); */
                        break;
                case 'I':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = s->seq[sp];
                                genome[rp] =255;
                                sp++;
                                rp++;
                        }
                        break;
                case 'D':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = 255;
                                rp++;
                        }
                        break;
                default:
                        break;
                }
        }
        aln_len = rp;

        /* for(int i = 0; i < aln_len;i++){ */
        /*         LOG_MSG("%d %d",genome[i],read[i]); */
        /* } */
        int l = a->md->len;
        if(l){
                int pos = 0;
                int i = 0;
                char tmp_num[8];
                char* md = a->md->str;
                while(i < l){
                        if(isdigit((int)md[i])){
                                int j = 0;
                                while (isdigit(md[i])) {
                                        tmp_num[j] = md[i];
                                        j++;
                                        i++;
                                        if(i == l){
                                                break;
                                        }
                                }
                                tmp_num[j] = 0;

                                int c = atoi(tmp_num);

                                //fprintf(stderr,"MD:%d\n",c);
                                for(j = 0; j < c;j++){
                                        if(genome[pos] != 255){
                                                genome[pos] = read[pos];
                                        }else{
                                                c++;
                                        }
                                        pos++;
                                }
                        }else if(isalpha((int)md[i])){
                                genome[pos] = (int)md[i];
                                pos++;
                                i++;
                        }else if(md[i] == '^'){
                                i++;
                        }
                        while(genome[pos] == 255){
                                pos++;
                        }
                }
        }
        /* LOG_MSG("GFINADSA : -------------------"); */
        /* for(int i = 0; i < aln_len;i++){ */
        /*         LOG_MSG("%d %d",genome[i],read[i]); */
        /* } */
        for(int i = 0; i < aln_len;i++){
                if(genome[i] == 255){
                        genome[i] = '-';
                }
                if(read[i] == 255){
                        read[i] = '-';
                }
        }
        char* cig = NULL;
        get_readable_cigar(a, &cig);
        fprintf(stdout,"Cigar: %s\n",cig);
        fprintf(stdout,"   MD: %s\n",TLD_STR(a->md));
        gfree(cig);

        fprintf(stdout, "%s (genome)\n",genome);
        for(int i = 0; i < aln_len;i++){
                if(genome[i] == read[i]){
                        fprintf(stdout,"|");
                }else{
                        fprintf(stdout," ");
                }
        }
        fprintf(stdout,"\n");

        fprintf(stdout, "%s (read)\n",read);
        gfree(genome);
        gfree(read);

        return OK;
ERROR:
        gfree(genome);
        gfree(read);
        return FAIL;
}

