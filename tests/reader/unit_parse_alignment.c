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

/* Per-record CIGAR counts (identity inputs, standard-plots.md 3.4) and the
 * mismatch marking of X ops when there is no MD (3.5). */
struct count_case {
        const char* name;
        const char* seq;
        const char* cigar;
        const char* md;         /* "" for none */
        int flag;
        uint32_t m, i_ev, i_bases, d_ev, d_bases, s_bases, x_ops;
        int has_md_or_x;
        const char* genome;     /* expected genome string, NULL: not checked */
};

static const struct count_case count_cases[] = {
        /* name         seq              cigar            md       flag m  iev ib dev db  s  x  has genome */
        /* one insertion (2 bases) and one deletion (1 base) */
        {"I_and_D",     "ACGTGGACGTAC",  "4M2I4M1D2M",    "8^T2",   0, 10, 1, 2, 1, 1, 0, 0, 1, NULL},
        /* two events of the same kind count as two events */
        {"two_I",       "ACGTGACGTAC",   "2M1I3M1I4M",    "9",      0, 9,  2, 2, 0, 0, 0, 0, 1, NULL},
        {"two_D",       "ACGTAC",        "2M3D2M2D2M",    "2^TTT2^GG2",0,6, 0, 0, 2, 5, 0, 0, 1, NULL},
        /* soft clips counted, hard clips and N / P ignored */
        {"soft_clip",   "TTTACGTACGGG",  "3S7M2S",        "7",      0, 7,  0, 0, 0, 0, 5, 0, 1, NULL},
        {"hard_clip",   "TTTACGTACG",    "5H3S7M4H",      "7",      0, 7,  0, 0, 0, 0, 3, 0, 1, NULL},
        {"N_and_P",     "ACGTACGTAC",    "5M100N2P5M",    "",       0, 10, 0, 0, 0, 0, 0, 0, 0, "ACGTACGTAC"},
        /* M only, no MD, no X: mismatches unknown */
        {"M_no_md",     "ACGTACGTAC",    "10M",           "",       0, 10, 0, 0, 0, 0, 0, 0, 0, "ACGTACGTAC"},
        /* X ops mark mismatches without MD; the reference base is unknown ('N') */
        {"X_no_md",     "ACGTACGTAC",    "5=1X4=",        "",       0, 10, 0, 0, 0, 0, 0, 1, 1, "ACGTANGTAC"},
        {"X_run_no_md", "ACGTACGTAC",    "2=3X5=",        "",       0, 10, 0, 0, 0, 0, 0, 3, 1, "ACNNNCGTAC"},
        /* a read base N in an X op still has to differ from the reference */
        {"X_readN_no_md","ACGTANGTAC",   "5=1X4=",        "",       0, 10, 0, 0, 0, 0, 0, 1, 1, "ACGTAAGTAC"},
        /* X ops in the alignment of a record with I and S */
        {"X_I_S_no_md", "TTACGTAGGAC",   "2S3=1X1I3=1S",  "",       0, 7,  1, 1, 0, 0, 3, 1, 1, "ACGN-GGA"},
        /* MD wins when present */
        {"X_with_md",   "ACGTACGTAC",    "5=1X4=",        "5G4",    0, 10, 0, 0, 0, 0, 0, 1, 1, "ACGTAGGTAC"},
        {"M_with_md",   "ACGTACGTAC",    "10M",           "4C5",    0, 10, 0, 0, 0, 0, 0, 0, 1, "ACGTCCGTAC"},
        /* unmapped: everything stays zero even after a record with counts */
        {"unmapped",    "ACGTACGTAC",    "10M",           "10",     4, 0,  0, 0, 0, 0, 0, 0, 0, NULL},
        /* CIGAR that does not fit SEQ: counts are not trusted */
        {"cigar_seq_mismatch","ACGTACGTAC","30M",         "30",     0, 0,  0, 0, 0, 0, 0, 0, 0, NULL},
        /* SEQ '*': nothing to align, but the CIGAR counts are still valid */
        {"seq_star",    "",              "10S40M",        "",       0, 40, 0, 0, 0, 0, 10, 0, 0, NULL},
        {"seq_star_X",  "",              "10=2X8=",       "",       0, 20, 0, 0, 0, 0, 0, 2, 1, NULL},
};

static int run_count_case(struct tl_seq_buffer* sb, const struct count_case* c)
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
        a->reverse = 0;
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
        CHECK(a->m_bases == c->m, "%s: m_bases %u, expected %u", c->name, a->m_bases, c->m);
        CHECK(a->i_ev == c->i_ev, "%s: i_ev %u, expected %u", c->name, a->i_ev, c->i_ev);
        CHECK(a->i_bases == c->i_bases, "%s: i_bases %u, expected %u", c->name, a->i_bases, c->i_bases);
        CHECK(a->d_ev == c->d_ev, "%s: d_ev %u, expected %u", c->name, a->d_ev, c->d_ev);
        CHECK(a->d_bases == c->d_bases, "%s: d_bases %u, expected %u", c->name, a->d_bases, c->d_bases);
        CHECK(a->s_bases == c->s_bases, "%s: s_bases %u, expected %u", c->name, a->s_bases, c->s_bases);
        CHECK(a->x_ops == c->x_ops, "%s: x_ops %u, expected %u", c->name, a->x_ops, c->x_ops);
        CHECK(a->has_md_or_x == c->has_md_or_x, "%s: has_md_or_x %d, expected %d", c->name, a->has_md_or_x, c->has_md_or_x);
        if(c->genome){
                size_t gl = strlen(c->genome);
                CHECK((size_t)a->aln_len == gl, "%s: aln_len %d, expected %zu", c->name, a->aln_len, gl);
                if((size_t)a->aln_len == gl){
                        CHECK(memcmp(a->genome, c->genome, gl) == 0,
                              "%s: genome %.*s, expected %s", c->name, a->aln_len, a->genome, c->genome);
                }
        }
        return OK;
}

/* a record without NM must report error == -1, also in a slot that held an NM before */
#define NM_FILE "unit_parse_alignment_nm.sam"
static int test_nm_reset(void)
{
        struct sam_bam_file* f = NULL;
        struct tl_seq_buffer* sb = NULL;
        FILE* fp = NULL;
        int ok = 0;

        fp = fopen(NM_FILE, "w");
        if(!fp){
                fprintf(stderr, "FAIL nm_reset: cannot write %s\n", NM_FILE);
                n_failed++;
                return OK;
        }
        fprintf(fp, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
        /* chunk 1 (4 slots): NM 3, none, NM 0, none */
        fprintf(fp, "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:3\n");
        fprintf(fp, "r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");
        fprintf(fp, "r3\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n");
        fprintf(fp, "r4\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");
        /* chunk 2 (same slots): none, NM 2, unmapped, none */
        fprintf(fp, "r5\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");
        fprintf(fp, "r6\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:2\n");
        fprintf(fp, "r7\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\t*\n");
        fprintf(fp, "r8\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");
        ok = (fclose(fp) == 0);
        if(!ok){
                fprintf(stderr, "FAIL nm_reset: cannot write %s\n", NM_FILE);
                n_failed++;
                remove(NM_FILE);
                return OK;
        }

        RUN(alloc_tl_seq_buffer(&sb, 4));
        RUN(add_aln_data(sb));
        RUN(open_bam(&f, NM_FILE));

        RUN(read_bam_chunk(f, sb));
        CHECK(sb->num_seq == 4, "nm_reset: chunk 1 has %d records", sb->num_seq);
        if(sb->num_seq == 4){
                const int expect1[4] = {3, -1, 0, -1};
                for(int i = 0; i < 4; i++){
                        struct aln_data* a = sb->sequences[i]->data;
                        CHECK(a->error == expect1[i], "nm_reset: chunk 1 record %d error %d, expected %d", i, a->error, expect1[i]);
                }
        }
        reset_tl_seq_buffer(sb);
        clear_aln_data(sb);
        RUN(read_bam_chunk(f, sb));
        CHECK(sb->num_seq == 4, "nm_reset: chunk 2 has %d records", sb->num_seq);
        if(sb->num_seq == 4){
                const int expect2[4] = {-1, 2, -1, -1};
                for(int i = 0; i < 4; i++){
                        struct aln_data* a = sb->sequences[i]->data;
                        CHECK(a->error == expect2[i], "nm_reset: chunk 2 record %d error %d, expected %d", i, a->error, expect2[i]);
                }
        }
        n_checks++;     /* the whole sequence ran */

        RUN(close_bam(f));
        f = NULL;
        RUN(free_aln_data(sb));
        free_tl_seq_buffer(sb);
        remove(NM_FILE);
        return OK;
ERROR:
        remove(NM_FILE);
        return FAIL;
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
        int nc = (int)(sizeof(count_cases) / sizeof(count_cases[0]));
        for(int i = 0; i < nc; i++){
                RUN(run_count_case(sb, &count_cases[i]));
        }
        RUN(free_aln_data(sb));
        free_tl_seq_buffer(sb);
        sb = NULL;
        RUN(test_nm_reset());
        n += nc + 1;
        fprintf(stdout, "unit_parse_alignment: %d cases, %d checks, %d failed\n", n, n_checks, n_failed);
        return n_failed ? 1 : 0;
ERROR:
        fprintf(stderr, "unit_parse_alignment: setup error\n");
        return 2;
}
