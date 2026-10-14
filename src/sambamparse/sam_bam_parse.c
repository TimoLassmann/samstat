
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
        a->n_clip5 = 0;
        a->n_clip3 = 0;
        a->m_bases = 0;
        a->i_ev = 0;
        a->i_bases = 0;
        a->d_ev = 0;
        a->d_bases = 0;
        a->s_bases = 0;
        a->x_ops = 0;
        a->has_md_or_x = 0;
        if(a->flag & BAM_FUNMAP){
                a->aln_len = 0;
                return OK;
        }
        aln_len = 0;
        int q_len = 0;  /* number of query bases consumed by the CIGAR (M, I, S, =, X) */

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
                        q_len += bam_cigar_oplen(a->cigar[i]);
                        a->s_bases += bam_cigar_oplen(a->cigar[i]);
                        break;
                case 'X':
                        a->x_ops += bam_cigar_oplen(a->cigar[i]);
                        /* fall through: an X op is also aligned bases */
                case 'M':
                case '=':
                        a->m_bases += bam_cigar_oplen(a->cigar[i]);
                        q_len += bam_cigar_oplen(a->cigar[i]);
                        aln_len += bam_cigar_oplen(a->cigar[i]);
                        break;
                case 'I':
                        a->i_ev++;
                        a->i_bases += bam_cigar_oplen(a->cigar[i]);
                        q_len += bam_cigar_oplen(a->cigar[i]);
                        aln_len += bam_cigar_oplen(a->cigar[i]);
                        break;
                case 'D':
                        Oplen = bam_cigar_oplen(a->cigar[i]);
                        a->d_ev++;
                        a->d_bases += (uint32_t)Oplen;
                        aln_len += Oplen;
                        break;
                default:
                        break;
                }
        }
        /* SEQ '*' (no bases) or a CIGAR that does not fit SEQ: nothing to align. */
        if(s->len == 0 || q_len != s->len){
                a->aln_len = 0;
                if(s->len != 0){
                        /* CIGAR disagrees with SEQ: the counts cannot be trusted */
                        a->m_bases = 0;
                        a->i_ev = 0;
                        a->i_bases = 0;
                        a->d_ev = 0;
                        a->d_bases = 0;
                        a->s_bases = 0;
                        a->x_ops = 0;
                }
                /* SEQ '*': the CIGAR counts stay valid, there is just nothing to align */
                a->has_md_or_x = (s->len == 0 && (a->x_ops > 0 || (a->md != NULL && a->md->len > 0))) ? 1 : 0;
                return OK;
        }

        /* A soft clip is the first (last) operation, ignoring hard clips. */
        int first = 0;
        int last = a->n_cigar - 1;
        while(first <= last && bam_cigar_op(a->cigar[first]) == BAM_CHARD_CLIP){
                first++;
        }
        while(last >= first && bam_cigar_op(a->cigar[last]) == BAM_CHARD_CLIP){
                last--;
        }
        if(first <= last){
                if(bam_cigar_op(a->cigar[first]) == BAM_CSOFT_CLIP){
                        a->n_clip5 = bam_cigar_oplen(a->cigar[first]);
                }
                if(bam_cigar_op(a->cigar[last]) == BAM_CSOFT_CLIP && last != first){
                        a->n_clip3 = bam_cigar_oplen(a->cigar[last]);
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
        rp = 0;
        sp = 0;
        int l = 0;
        if(a->md != NULL){
                l = a->md->len;
        }
        /* MD wins; without it X ops are the only source of mismatches */
        a->has_md_or_x = (l > 0 || a->x_ops > 0) ? 1 : 0;

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
                        sp += Oplen;
                        break;
                case 'M':
                case '=':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = seq[sp];
                                genome[rp] = read[rp];
                                sp++;
                                rp++;
                        }
                        break;
                case 'X':
                        for(int j = 0; j < Oplen;j++){
                                read[rp] = seq[sp];
                                genome[rp] = read[rp];
                                if(!l){
                                        /* no MD: the reference base is unknown but differs */
                                        genome[rp] = (read[rp] == 'N') ? 'A' : 'N';
                                }
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
        if(l){
                int pos = 0;
                int i = 0;
                char* md = TLD_STR(a->md); //->str;
                /* Everything below is bounded by the alignment length: an MD tag
                 * that disagrees with the CIGAR must not reach outside genome / read. */
                while(i < l && pos < aln_len){
                        if(isdigit((int)md[i])){
                                int c = 0;
                                while(i < l && isdigit((int)md[i])){
                                        /* saturate - anything above aln_len covers the whole alignment */
                                        if(c <= aln_len){
                                                c = c * 10 + (md[i] - '0');
                                        }
                                        i++;
                                }
                                /* c matches; insertions are not part of the reference */
                                while(c && pos < aln_len){
                                        if(genome[pos] != '-'){
                                                genome[pos] = read[pos];
                                                c--;
                                        }
                                        pos++;
                                }
                        }else if(isalpha((int)md[i])){
                                genome[pos] = (int)md[i];
                                pos++;
                                i++;
                        }else{
                                /* '^' or anything unexpected */
                                i++;
                        }
                        while(pos < aln_len && genome[pos] == '-'){
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
