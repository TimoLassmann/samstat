
#include "sam.h"
#include "tld.h"
#include "../htsinterface/htsglue.h"
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
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
        uint8_t* seq = NULL;
        int sp;
        int rp;
        char Op;
        int Oplen;
        int aln_len;

        a = s->data;
        if(a->flag & BAM_FUNMAP){
                a->aln_len = 0;
                return OK;
        }
        aln_len = 0;

        for(int i = 0;i < a->n_cigar;i++){
                Op = bam_cigar_opchr(a->cigar[i]);
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
                        if(i == 0){
                                a->n_clip5 = bam_cigar_oplen(a->cigar[i]);
                        }
                        if(i == a->n_cigar - 1){
                                a->n_clip3 = bam_cigar_oplen(a->cigar[i]);
                        }
                        break;
                case 'M':
                case '=':
                case 'X':
                case 'I':
                case 'D':
                        Oplen = bam_cigar_oplen(a->cigar[i]);
                        aln_len += Oplen;
                        break;
                default:
                        break;
                }
        }
        a->aln_len = aln_len;
        aln_len++;

        if(a->aln_len >= a->aln_len_alloc){
                a->aln_len_alloc = MACRO_MAX(a->aln_len_alloc + a->aln_len_alloc / 2, a->aln_len_alloc + a->aln_len - a->aln_len_alloc + 1);
                gfree(a->read);
                gfree(a->genome);
                a->read = NULL;
                a->genome = NULL;
                galloc(&a->genome, a->aln_len_alloc);
                galloc(&a->read, a->aln_len_alloc);
        }

        genome = a->genome;
        read = a->read;
        /* for(int i = 0;i < aln_len;i++){ */
        /*         genome[i] = '-'; */
        /*         read[i] = '-'; */
        /* } */
        rp = 0;
        sp = 0;

        seq = s->seq->str;
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
                                read[rp] = seq[sp];
                                genome[rp] = read[rp];
                                sp++;
                                rp++;
                        }
                        break;
                case 'I':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = seq[sp];
                                genome[rp] = '-';
                                sp++;
                                rp++;
                        }
                        break;
                case 'D':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = '-';
                                genome[rp] = 'N';
                                rp++;
                        }
                        break;
                default:
                        break;
                }
        }
        /* LOG_MSG("RP: %d", rp); */
        aln_len = rp;
        /* fprintf(stdout,"%s - genome\n%s - read\n", genome , read); */
        int l = 0;
        if(a->md == NULL){
                l = 0;
        }else{
                l = a->md->len;
        }
        if(l){
                int pos = 0;
                int i = 0;
                char tmp_num[8];
                char* md = TLD_STR(a->md); //->str;
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
                                for(j = 0; j < c;j++){
                                        if(genome[pos] != '-'){
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
                        while(genome[pos] == '-'){
                                pos++;
                        }
                }
        }/* else{ */
        /*         for(int i = 0; i < a->aln_len;i++){ */
        /*                 if(read[i] == 255){ */
        /*                         genome[i] = 'N'; */
        /*                 }else{ */
        /*                         genome[i] = read[i]; */
        /*                 } */
        /*                 /\* fprintf(stdout,"%c",genome[i]); *\/ */
        /*         } */
        /*         fprintf(stdout,"Filing done \n"); */
        /* } */

        /* for(int i = 0; i < aln_len;i++){ */
        /*         if(genome[i] == 255){ */
        /*                 genome[i] = '-'; */
        /*         } */
        /*         if(read[i] == 255){ */
        /*                 read[i] = '-'; */
        /*         } */
        /* } */
        /* fprintf(stdout,"%s - genome len: %d\n%s - read\n", genome ,aln_len, read); */
        return OK;
ERROR:
        return FAIL;
}

int fix_orientation(struct tl_seq *s)
{
        struct aln_data* a = NULL;

        a = s->data;
        if(a->reverse){
                /* LOG_MSG("Oh no I should reverse!"); */
                rev_comp_tl_seq(s);
                if(a->aln_len){
                        reverse_comp(a->genome, a->aln_len);
                        reverse_comp(a->read, a->aln_len);
                }
                a->reverse = 0;
                uint32_t tmp;
                tmp = a->n_clip5;
                a->n_clip5 = a->n_clip3;
                a->n_clip3 = tmp;
        }
        return OK;
}
