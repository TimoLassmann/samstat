
#include "sam.h"
#include "tldlite.h"
#include "../htsinterface/htsglue.h"
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define SAM_BAM_PARSE_IMPORT
#include "sam_bam_parse.h"


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
        int64_t aln_len64;
        uint64_t d_sum;
        uint64_t m_sum;
        uint64_t i_sum;
        uint64_t s_sum;
        uint64_t x_sum;

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
        /*
         * 64-bit sums: a CIGAR operation can be 268 million bases long and a record can have many.
         * The alignment strings get one column per read base in M, =, X and I, and ONE column per
         * D operation, however long: their length is bounded by the read length plus the number of
         * operations, never by the reference span (a 100 Mb deletion must not need 200 Mb).
         */
        aln_len64 = 0;
        d_sum = 0;
        m_sum = 0;
        i_sum = 0;
        s_sum = 0;
        x_sum = 0;
        int64_t q_len = 0;  /* number of query bases consumed by the CIGAR (M, I, S, =, X) */
        int eqx = 0;        /* the CIGAR uses = or X: every aligned base is a known match or mismatch */

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
                        s_sum += bam_cigar_oplen(a->cigar[i]);
                        break;
                case 'X':
                        x_sum += bam_cigar_oplen(a->cigar[i]);
                        /* fall through: an X op is also aligned bases */
                case 'M':
                case '=':
                        if(Op != 'M'){
                                eqx = 1;
                        }
                        m_sum += bam_cigar_oplen(a->cigar[i]);
                        q_len += bam_cigar_oplen(a->cigar[i]);
                        aln_len64 += bam_cigar_oplen(a->cigar[i]);
                        break;
                case 'I':
                        a->i_ev++;
                        i_sum += bam_cigar_oplen(a->cigar[i]);
                        q_len += bam_cigar_oplen(a->cigar[i]);
                        aln_len64 += bam_cigar_oplen(a->cigar[i]);
                        break;
                case 'D':
                        a->d_ev++;
                        d_sum += bam_cigar_oplen(a->cigar[i]);
                        aln_len64 += 1;
                        break;
                default:
                        break;
                }
        }
        /* saturated: a larger sum makes the NM-based identity negative, i.e. not available.
         * The other sums are bounded by the read length when the CIGAR fits SEQ, but a
         * SEQ '*' record keeps its counts and its 28-bit operations can add up past 2^32. */
        a->d_bases = (d_sum > UINT32_MAX) ? UINT32_MAX : (uint32_t) d_sum;
        a->m_bases = (m_sum > UINT32_MAX) ? UINT32_MAX : (uint32_t) m_sum;
        a->i_bases = (i_sum > UINT32_MAX) ? UINT32_MAX : (uint32_t) i_sum;
        a->s_bases = (s_sum > UINT32_MAX) ? UINT32_MAX : (uint32_t) s_sum;
        a->x_ops = (x_sum > UINT32_MAX) ? UINT32_MAX : (uint32_t) x_sum;
        if(q_len != (int64_t) s->len || aln_len64 > INT32_MAX){
                /* treated below as a CIGAR that does not fit SEQ */
                q_len = -1;
                aln_len64 = 0;
        }
        aln_len = (int) aln_len64;
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
                a->has_md_or_x = (s->len == 0 && (eqx || (a->md != NULL && a->md->len > 0))) ? 1 : 0;
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
        /* set only once the buffers hold aln_len + 1 bytes: after a failed allocation the
         * record has nothing aligned (aln_len 0) and no buffer (aln_len_alloc 0), and the
         * next record allocates again */
        a->aln_len = 0;
        if(aln_len >= a->aln_len_alloc){
                const int64_t grow = (int64_t) a->aln_len_alloc + a->aln_len_alloc / 2;
                const int64_t need = (int64_t) aln_len + 1;
                int64_t new_alloc = MACRO_MAX(grow, need);

                ASSERT(need <= INT32_MAX, "Alignment of %s is too long", TLD_STR(s->name));
                if(new_alloc > INT32_MAX){
                        new_alloc = INT32_MAX;
                }
                if(a->read){
                        gfree(a->read);
                }
                if(a->genome){
                        gfree(a->genome);
                }
                a->read = NULL;
                a->genome = NULL;
                a->aln_len_alloc = 0;
                galloc(&a->genome, (int) new_alloc);
                galloc(&a->read, (int) new_alloc);
                a->aln_len_alloc = (int) new_alloc;
        }
        a->aln_len = aln_len;
        aln_len++;

        genome = a->genome;
        read = a->read;
        rp = 0;
        sp = 0;
        int l = 0;
        if(a->md != NULL){
                l = a->md->len;
        }
        /* MD wins; without it a CIGAR with = and X ops says which bases match. A read of = ops only
         * matches perfectly: its mismatches are known (none), so it counts in the mismatch rate */
        a->has_md_or_x = (l > 0 || eqx) ? 1 : 0;

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
                        /* one column per deletion, whatever its length (see above) */
                        read[rp] = '-';
                        genome[rp] = 'N';
                        rp++;
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
                        if(isdigit((unsigned char)md[i])){
                                /* 64 bit and saturated: c <= aln_len <= INT32_MAX before the
                                 * multiply, so c * 10 + 9 cannot overflow */
                                int64_t c = 0;
                                while(i < l && isdigit((unsigned char)md[i])){
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
                        }else if(md[i] == '^'){
                                /* deleted reference bases: all of them belong to the one column of
                                 * the deletion, which shows the first of them */
                                int first_del = 0;
                                i++;
                                while(i < l && isalpha((unsigned char)md[i])){
                                        if(!first_del){
                                                first_del = md[i];
                                        }
                                        i++;
                                }
                                if(read[pos] == '-'){
                                        if(first_del){
                                                genome[pos] = first_del;
                                        }
                                        pos++;
                                }
                        }else if(isalpha((unsigned char)md[i])){
                                genome[pos] = (uint8_t)md[i];
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
