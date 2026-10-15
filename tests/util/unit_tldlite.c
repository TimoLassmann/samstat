/* Unit test for tldlite.h: string buffers, galloc, file names, reverse complement and
 * sequence buffers. Small inputs only; the growth check appends 100 000 characters. */
#include "tldlite.h"

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

static int test_strbuf(void)
{
        tld_strbuf* b = NULL;
        RUN(tld_strbuf_alloc(&b, 4));
        CHECK(b->len == 0 && b->alloc_len == 4 && TLD_STR(b)[0] == 0, "new buffer");
        RUN(tld_append(b, "abc"));
        RUN(tld_append_char(b, 'd'));
        RUN(tld_append(b, "efghij"));
        CHECK(strcmp(TLD_STR(b), "abcdefghij") == 0 && TLD_STRLEN(b) == 10, "appends: %s", TLD_STR(b));
        CHECK(b->alloc_len > b->len, "room for the terminating 0 (%d of %d)", b->len, b->alloc_len);
        RUN(tld_strbuf_clear(b));
        CHECK(b->len == 0 && TLD_STR(b)[0] == 0, "clear");

        /* growth doubles the capacity: n appends cost O(n), and the capacity records the
         * real size (the tld bug fixed upstream in c32f857 left alloc_len unchanged) */
        int grows = 0;
        int cap = b->alloc_len;
        for(int i = 0; i < 100000; i++){
                RUN(tld_append_char(b, (char) ('a' + i % 26)));
                if(b->alloc_len != cap){
                        grows++;
                        cap = b->alloc_len;
                }
        }
        CHECK(b->len == 100000 && b->alloc_len > b->len, "100 000 appends: len %d alloc %d", b->len, b->alloc_len);
        CHECK(grows <= 20, "capacity grew %d times", grows);
        CHECK(TLD_STR(b)[0] == 'a' && TLD_STR(b)[99999] == (char) ('a' + 99999 % 26) && TLD_STR(b)[100000] == 0, "content kept");
        tld_strbuf_free(b);
        b = NULL;

        RUN(tld_strbuf_alloc(&b, 0));
        CHECK(b->alloc_len == 1024, "size 0 gives 1024 (%d)", b->alloc_len);
        tld_strbuf_free(b);
        tld_strbuf_free(NULL);
        return OK;
ERROR:
        tld_strbuf_free(b);
        return FAIL;
}

static int test_galloc(void)
{
        uint64_t* v = NULL;
        uint64_t** m = NULL;
        char** names = NULL;

        galloc(&v, 7);
        int zero = 1;
        for(int i = 0; i < 7; i++){
                zero &= (v[i] == 0);
                v[i] = (uint64_t) i;
        }
        CHECK(zero && v[6] == 6, "1D: zeroed and writable");

        galloc(&m, 3, 5);
        zero = 1;
        for(int i = 0; i < 3; i++){
                for(int j = 0; j < 5; j++){
                        zero &= (m[i][j] == 0);
                        m[i][j] = (uint64_t) (10 * i + j);
                }
        }
        CHECK(zero, "2D: zeroed");
        CHECK(m[1] == m[0] + 5 && m[2] == m[0] + 10, "2D: rows in one block");
        CHECK(m[2][4] == 24 && m[0][9] == 14, "2D: row-major values");

        galloc(&names, 2, 16);
        snprintf(names[0], 16, "first");
        snprintf(names[1], 16, "second");
        CHECK(strcmp(names[0], "first") == 0 && strcmp(names[1], "second") == 0, "2D char");

        gfree(v);
        gfree(m);
        gfree(names);
        return OK;
ERROR:
        return FAIL;
}

static int test_files(void)
{
        char* f = NULL;
        RUN(tld_get_filename("dir/sub/reads.bam", &f));
        CHECK(strcmp(f, "reads.bam") == 0, "path: %s", f);
        MFREE(f);
        RUN(tld_get_filename("reads.fq.gz", &f));
        CHECK(strcmp(f, "reads.fq.gz") == 0, "no directory: %s", f);
        MFREE(f);
        printf("(an ERROR line for \"dir/\" is expected)\n");
        CHECK(tld_get_filename("dir/", &f) == FAIL && f == NULL, "no file name");

        CHECK(tld_file_exists(".") == OK, "existing path");
        CHECK(tld_file_exists("/no/such/samstat/file") == FAIL, "missing path");
        return OK;
ERROR:
        return FAIL;
}

static int test_revcomp(void)
{
        uint8_t odd[] = "ACGTN";
        uint8_t even[] = "aacgtT";
        uint8_t one[] = "G";
        uint8_t q[] = "12345";

        reverse_comp(odd, 5);
        CHECK(memcmp(odd, "NACGT", 5) == 0, "odd: %s", odd);
        reverse_comp(even, 6);
        CHECK(memcmp(even, "Aacgtt", 6) == 0, "even, lower case: %s", even);
        reverse_comp(one, 1);
        CHECK(one[0] == 'C', "one base: %c", one[0]);
        reverse_comp(one, 0);
        CHECK(one[0] == 'C', "zero length leaves the buffer");
        reverse(q, 5);
        CHECK(memcmp(q, "54321", 5) == 0, "reverse: %s", q);
        return OK;
}

static int test_seq_buffer(void)
{
        struct tl_seq_buffer* sb = NULL;

        RUN(alloc_tl_seq_buffer(&sb, 3));
        CHECK(sb->malloc_num == 3 && sb->num_seq == 0 && sb->data == NULL, "new buffer");
        CHECK(((uintptr_t) sb) % 64 == 0 && ((uintptr_t) sb->sequences[2]) % 64 == 0, "64-byte aligned");

        struct tl_seq* s = sb->sequences[0];
        RUN(tld_append(s->seq, "AACG"));
        RUN(tld_append(s->qual, "ABCD"));
        RUN(tld_append(s->name, "r1"));
        s->len = 4;
        sb->num_seq = 1;
        RUN(rev_comp_tl_seq(s));
        CHECK(strcmp(TLD_STR(s->seq), "CGTT") == 0 && strcmp(TLD_STR(s->qual), "DCBA") == 0,
              "rev_comp_tl_seq: %s %s", TLD_STR(s->seq), TLD_STR(s->qual));

        RUN(reset_tl_seq_buffer(sb));
        CHECK(sb->num_seq == 0 && s->len == 0 && s->seq->len == 0 && s->qual->len == 0 && s->name->len == 0, "reset");
        free_tl_seq_buffer(sb);
        free_tl_seq_buffer(NULL);
        return OK;
ERROR:
        free_tl_seq_buffer(sb);
        return FAIL;
}

int main(void)
{
        CHECK(test_strbuf() == OK, "string buffers");
        CHECK(test_galloc() == OK, "galloc");
        CHECK(test_files() == OK, "files");
        CHECK(test_revcomp() == OK, "reverse complement");
        CHECK(test_seq_buffer() == OK, "sequence buffers");
        printf("unit_tldlite: %d checks, %d failed\n", n_checks, n_failed);
        return n_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
