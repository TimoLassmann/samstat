/* Unit test for parse_alignment() / fix_orientation().
 *
 * A struct tl_seq and its aln_data are filled by hand (SEQ, CIGAR, MD, strand)
 * and the produced read / genome strings and clip lengths are compared with
 * values worked out from the SAM specification.
 */
#include "tld.h"
#include "seq/tld-seq.h"
#include "../../src/htsinterface/htsglue.h"
#include "../../src/sambamparse/sam_bam_parse.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n_checks = 0;
static int n_failed = 0;

#define CHECK(COND, ...) do{                                            \
                n_checks++;                                             \
                if(!(COND)){                                            \
                        n_failed++;                                     \
                        fprintf(stderr, "FAIL (%s:%d): %s: ", __FILE__, __LINE__, #COND); \
                        fprintf(stderr, __VA_ARGS__);                   \
                        fprintf(stderr, "\n");                          \
                }                                                       \
        }while(0)

struct aln_case {
        const char* name;
        const char* seq;        /* SEQ as stored in the record ("" for '*') */
        const char* cigar;      /* e.g. "5H3S7M" */
        const char* md;         /* MD tag value, "" for none */
        int reverse;
        int flag;
        /* expected values after parse_alignment() */
        int aln_len;
        const char* read;
        const char* genome;
        int clip5;
        int clip3;
        /* expected values after fix_orientation() (reverse strand only) */
        const char* rc_seq;
};

static int set_cigar(struct aln_data* a, const char* cigar)
{
        const char* ops = BAM_CIGAR_STR;
        int n = 0;
        const char* p = cigar;
        while(*p){
                char* end = NULL;
                long len = strtol(p, &end, 10);
                const char* op = strchr(ops, *end);
                if(end == p || !op){
                        return FAIL;
                }
                if(n >= a->n_alloc_cigar){
                        return FAIL;
                }
                a->cigar[n++] = bam_cigar_gen((uint32_t)len, (uint32_t)(op - ops));
                p = end + 1;
        }
        a->n_cigar = n;
        return OK;
}

static int run_case(struct tl_seq_buffer* sb, const struct aln_case* c)
{
        struct tl_seq* s = NULL;
        struct aln_data* a = NULL;

        reset_tl_seq_buffer(sb);
        clear_aln_data(sb);
        s = sb->sequences[0];
        a = s->data;

        tld_append(s->seq, (char*)c->seq);
        tld_append(s->qual, "");
        s->len = (int)strlen(c->seq);
        if(c->md[0]){
                tld_append(a->md, (char*)c->md);
        }
        a->flag = (uint16_t)c->flag;
        a->reverse = (uint8_t)c->reverse;
        if(set_cigar(a, c->cigar) != OK){
                fprintf(stderr, "FAIL %s: bad CIGAR %s\n", c->name, c->cigar);
                n_failed++;
                return OK;
        }

        if(parse_alignment(s) != OK){
                fprintf(stderr, "FAIL %s: parse_alignment returned FAIL\n", c->name);
                n_failed++;
                return OK;
        }
        CHECK(a->aln_len == c->aln_len, "%s: aln_len %d, expected %d", c->name, a->aln_len, c->aln_len);
        CHECK(a->n_clip5 == (uint32_t)c->clip5, "%s: n_clip5 %u, expected %d", c->name, a->n_clip5, c->clip5);
        CHECK(a->n_clip3 == (uint32_t)c->clip3, "%s: n_clip3 %u, expected %d", c->name, a->n_clip3, c->clip3);
        if(a->aln_len == c->aln_len && a->aln_len > 0){
                CHECK(memcmp(a->read, c->read, (size_t)a->aln_len) == 0,
                      "%s: read %.*s, expected %s", c->name, a->aln_len, a->read, c->read);
                CHECK(memcmp(a->genome, c->genome, (size_t)a->aln_len) == 0,
                      "%s: genome %.*s, expected %s", c->name, a->aln_len, a->genome, c->genome);
        }

        if(c->reverse){
                RUN(fix_orientation(s));
                CHECK(a->reverse == 0, "%s: reverse flag not cleared", c->name);
                /* orientation fix swaps clips and reverse complements everything */
                CHECK(a->n_clip5 == (uint32_t)c->clip3, "%s: after fix n_clip5 %u, expected %d", c->name, a->n_clip5, c->clip3);
                CHECK(a->n_clip3 == (uint32_t)c->clip5, "%s: after fix n_clip3 %u, expected %d", c->name, a->n_clip3, c->clip5);
                CHECK(strcmp((char*)s->seq->str, c->rc_seq) == 0,
                      "%s: seq %s, expected %s", c->name, (char*)s->seq->str, c->rc_seq);
        }
        return OK;
ERROR:
        return FAIL;
}

/* Reverse complements checked by hand:
 *   read ACGTACG    -> CGTACGT     genome ACTTACG -> CGTAAGT
 *   seq TTTACGTACG  -> CGTACGTAAA                                        */
static const struct aln_case cases[] = {
        /* name        seq           cigar        md      rev flag  aln read         genome       c5 c3 rc_seq */
        {"M",          "ACGTACGTAC", "10M",       "10",    0, 0,   10, "ACGTACGTAC", "ACGTACGTAC", 0, 0, NULL},
        {"M_mismatch", "ACGTACGTAC", "10M",       "4C5",   0, 0,   10, "ACGTACGTAC", "ACGTCCGTAC", 0, 0, NULL},
        {"eq_x",       "ACGTACGTAC", "5=1X4=",    "5G4",   0, 0,   10, "ACGTACGTAC", "ACGTAGGTAC", 0, 0, NULL},
        {"I",          "ACGTGGACGT", "4M2I4M",    "8",     0, 0,   10, "ACGTGGACGT", "ACGT--ACGT", 0, 0, NULL},
        {"D",          "ACGTACGT",   "4M2D4M",    "4^TT4", 0, 0,   10, "ACGT--ACGT", "ACGTTTACGT", 0, 0, NULL},
        {"mismatch_after_I", "ACGTGGACGT", "4M2I4M", "6A1", 0, 0,  10, "ACGTGGACGT", "ACGT--ACAT", 0, 0, NULL},
        {"S5",         "TTTACGTACG", "3S7M",      "7",     0, 0,    7, "ACGTACG",    "ACGTACG",    3, 0, NULL},
        {"S5_mismatch","TTTACGTACG", "3S7M",      "2T4",   0, 0,    7, "ACGTACG",    "ACTTACG",    3, 0, NULL},
        {"S3",         "ACGTACGTTT", "7M3S",      "7",     0, 0,    7, "ACGTACG",    "ACGTACG",    0, 3, NULL},
        {"S5_S3",      "TTACGTACGGG","2S7M2S",    "2T4",   0, 0,    7, "ACGTACG",    "ACTTACG",    2, 2, NULL},
        {"H_S",        "TTTACGTACGTT","5H3S7M2S4H","7",    0, 0,    7, "ACGTACG",    "ACGTACG",    3, 2, NULL},
        {"H_only",     "ACGTACG",    "4H7M3H",    "7",     0, 0,    7, "ACGTACG",    "ACGTACG",    0, 0, NULL},
        {"S_after_H_5","TTTACGTACG", "5H3S7M",    "7",     0, 0,    7, "ACGTACG",    "ACGTACG",    3, 0, NULL},
        {"S_before_H_3","ACGTACGTTT","7M3S5H",    "7",     0, 0,    7, "ACGTACG",    "ACGTACG",    0, 3, NULL},
        {"rev_M",      "ACGTACG",    "7M",        "2T4",   1, 16,   7, "ACGTACG",    "ACTTACG",    0, 0, "CGTACGT"},
        {"rev_S5",     "TTTACGTACG", "3S7M",      "2T4",   1, 16,   7, "ACGTACG",    "ACTTACG",    3, 0, "CGTACGTAAA"},
        {"rev_H_S",    "TTTACGTACGTT","5H3S7M2S4H","7",    1, 16,   7, "ACGTACG",    "ACGTACG",    3, 2, "AACGTACGTAAA"},
        /* robustness: nothing to align, must not touch memory outside the buffers */
        {"seq_star",   "",           "10S40M",    "",      0, 0,    0, NULL,         NULL,         0, 0, NULL},
        {"seq_star_MD","",           "10S40M",    "40",    0, 0,    0, NULL,         NULL,         0, 0, NULL},
        {"unmapped",   "ACGTACGTAC", "10M",       "10",    0, 4,    0, NULL,         NULL,         0, 0, NULL},
        {"cigar_seq_mismatch", "ACGTACGTAC", "30M","30",   0, 0,    0, NULL,         NULL,         0, 0, NULL},
        /* MD that disagrees with the alignment */
        {"MD_too_long","ACGTACGTAC", "10M",       "1000",  0, 0,   10, "ACGTACGTAC", "ACGTACGTAC", 0, 0, NULL},
        {"MD_mm_then_long","ACGTACGTAC","10M",    "4C5000",0, 0,   10, "ACGTACGTAC", "ACGTCCGTAC", 0, 0, NULL},
        {"MD_8_digits","ACGTACGTAC", "10M",       "00000010",0, 0, 10, "ACGTACGTAC", "ACGTACGTAC", 0, 0, NULL},
        {"MD_30_digits","ACGTACGTAC","10M",       "999999999999999999999999999999",0, 0, 10, "ACGTACGTAC", "ACGTACGTAC", 0, 0, NULL},
        {"MD_letters_long","ACGTACGTAC","10M",    "ACGTACGTACGTACGTACGT",0, 0, 10, "ACGTACGTAC", "ACGTACGTAC", 0, 0, NULL},
        {"MD_del_long","ACGTACGT",   "4M2D4M",    "4^TTTTTTTT4",0, 0, 10, "ACGT--ACGT", "ACGTTTTTTT", 0, 0, NULL},
        {"MD_garbage", "ACGTACGTAC", "10M",       "4*#5",  0, 0,   10, "ACGTACGTAC", "ACGTACGTAC", 0, 0, NULL},
};

int main(void)
{
        struct tl_seq_buffer* sb = NULL;
        int n = (int)(sizeof(cases) / sizeof(cases[0]));

        RUN(alloc_tl_seq_buffer(&sb, 4));
        RUN(add_aln_data(sb));
        for(int i = 0; i < n; i++){
                RUN(run_case(sb, &cases[i]));
        }
        RUN(free_aln_data(sb));
        free_tl_seq_buffer(sb);

        fprintf(stdout, "unit_parse_alignment: %d cases, %d checks, %d failed\n", n, n_checks, n_failed);
        return n_failed ? 1 : 0;
ERROR:
        fprintf(stderr, "unit_parse_alignment: setup error\n");
        return 2;
}
