/* Unit test for tldlite.h: messages, macros, string buffers, galloc, file names, reverse
 * complement and sequence buffers. Small inputs only; the growth check appends 100 000
 * characters. The message checks write two small files in the working directory (the build
 * tree when run by ctest). */
#include "tldlite.h"

#include <unistd.h>

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
        v = NULL;
        m = NULL;
        names = NULL;

        /* the data block holds rows x cols elements of the element type, not of the pointer
         * type: on 32-bit, uint64_t elements overflowed a block sized with sizeof(pointer)
         * (the sanitizer build catches the write at the last index) */
        galloc(&m, 7, 9);
        m[6][8] = UINT64_MAX;
        m[0][0] = 1;
        CHECK(m[6][8] == UINT64_MAX && m[0][0] == 1 && m[6] == m[0] + 6 * 9, "2D uint64_t: last element");
        galloc(&names, 5, 3);
        names[4][2] = 'z';
        names[0][0] = 'a';
        CHECK(names[4][2] == 'z' && names[4] == names[0] + 12 && names[3][2] == 0, "2D char: last element");
        gfree(m);
        gfree(names);

        /* gfree / tld_free of NULL is a no-op */
        gfree(NULL);
        tld_free(NULL);
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
        fprintf(stderr, "(an error line for \"dir/\" is expected)\n");
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

/* --- messages and macros ------------------------------------------------ */

#define MSG_ERR_FILE "unit_tldlite_stderr.txt"
#define MSG_OUT_FILE "unit_tldlite_stdout.txt"

static int fail_silently(void)
{
        return FAIL;
}

static int succeed(void)
{
        return OK;
}

/* RUN of a failing function prints nothing (the root cause prints its own ERROR_MSG) */
static int run_fails(void)
{
        RUN(succeed());
        RUN(fail_silently());
        return OK;
ERROR:
        return FAIL;
}

/* ASSERT prints only its message: the test text (here with a '%') is not a format */
static int assert_fails(int n)
{
        ASSERT(n % 2 == 0, "odd value %d", n);
        return OK;
ERROR:
        return FAIL;
}

static int error_msg_fails(void)
{
        ERROR_MSG("bad value %d in %s", 7, "x.sam");
        return OK;
ERROR:
        return FAIL;
}

static int n_size_calls = 0;

static size_t next_size(void)
{
        n_size_calls++;
        return 32;
}

/* The macros are single statements: they need the ';' and work as an unbraced if / else body
 * (with the old if-without-else / trailing-';' forms this did not compile). The size argument
 * of MMALLOC is evaluated once. */
static int macro_statements(int flag)
{
        char* p = NULL;
        char* q = NULL;
        uint64_t* v = NULL;
        uint64_t** m = NULL;

        if(flag) RUN(succeed()); else RUN(fail_silently());
        if(flag) ASSERT(flag == 1, "flag %d", flag); else ERROR_MSG("not reached");
        if(flag) MMALLOC(p, next_size()); else MMALLOC(q, next_size());
        if(flag) galloc(&v, 4); else galloc(&m, 4, 4);
        CHECK(n_size_calls == 1, "MMALLOC evaluated its size %d times", n_size_calls);
        CHECK(p != NULL && p[31] == 0 && v != NULL && v[3] == 0, "allocated and zeroed");
        if(flag) MFREE(p); else MFREE(q);
        CHECK(p == NULL, "MFREE sets the pointer to NULL");
        /* MFREE of a NULL pointer is a silent no-op */
        MFREE(p);
        MFREE(q);
        gfree(v);
        gfree(m);
        return OK;
ERROR:
        MFREE(p);
        MFREE(q);
        gfree(v);
        gfree(m);
        return FAIL;
}

/* Strings that a near-full buffer cannot take fail cleanly instead of overflowing the int32_t
 * length. The buffer is faked (len and alloc_len near INT32_MAX over a small allocation):
 * the check runs before any memory is touched, so nothing large is allocated. */
static int strbuf_limit(void)
{
        tld_strbuf* b = NULL;
        int32_t len = 0;
        int32_t alloc_len = 0;
        int r = OK;

        RUN(tld_strbuf_alloc(&b, 16));
        len = b->len;
        alloc_len = b->alloc_len;
        b->alloc_len = INT32_MAX;
        b->len = INT32_MAX - 1;
        r = tld_append(b, "ab");
        CHECK(r == FAIL && b->len == INT32_MAX - 1 && b->alloc_len == INT32_MAX, "append past INT32_MAX fails");
        r = tld_append_char(b, 'a');
        CHECK(r == FAIL && b->len == INT32_MAX - 1, "append_char past INT32_MAX fails");
        b->len = len;
        b->alloc_len = alloc_len;
        tld_strbuf_free(b);
        return OK;
ERROR:
        tld_strbuf_free(b);
        return FAIL;
}

static char* read_file(const char* name, size_t* n)
{
        FILE* f = NULL;
        char* buf = NULL;
        long sz = 0;

        *n = 0;
        if((f = fopen(name, "rb")) == NULL){
                return NULL;
        }
        if(fseek(f, 0, SEEK_END) == 0 && (sz = ftell(f)) >= 0 && fseek(f, 0, SEEK_SET) == 0){
                buf = malloc((size_t) sz + 1);
                if(buf != NULL){
                        *n = fread(buf, 1, (size_t) sz, f);
                        buf[*n] = 0;
                }
        }
        fclose(f);
        return buf;
}

/* LOG / WARNING / ERROR go to stderr (never stdout), one line each, as
 * "samstat: [warning: |error: ]<msg>", without time stamps or source locations. stderr and
 * stdout are redirected to files while the messages are written, then restored. */
static int test_messages(void)
{
        char long_msg[3001];
        char* err = NULL;
        char* out = NULL;
        size_t n_err = 0;
        size_t n_out = 0;
        int saved_err = -1;
        int saved_out = -1;
        int r[6] = {0};

        memset(long_msg, 'x', 3000);
        long_msg[3000] = 0;

        fflush(stdout);
        fflush(stderr);
        saved_err = dup(fileno(stderr));
        saved_out = dup(fileno(stdout));
        if(saved_err < 0 || saved_out < 0 ||
           freopen(MSG_ERR_FILE, "w", stderr) == NULL || freopen(MSG_OUT_FILE, "w", stdout) == NULL){
                if(saved_err >= 0){
                        dup2(saved_err, fileno(stderr));
                }
                if(saved_out >= 0){
                        dup2(saved_out, fileno(stdout));
                }
                CHECK(0, "cannot redirect stderr / stdout");
                return FAIL;
        }

        LOG_MSG("hello %s", "world");
        WARNING_MSG("w %d", 1);
        r[0] = error_msg_fails();
        r[1] = assert_fails(3);
        r[2] = run_fails();
        r[3] = macro_statements(1);
        r[4] = strbuf_limit();
        WARNING_MSG("%s", long_msg);
        gfree(NULL);
        LOG_MSG("100%% done");

        fflush(stdout);
        fflush(stderr);
        dup2(saved_err, fileno(stderr));
        dup2(saved_out, fileno(stdout));
        close(saved_err);
        close(saved_out);
        clearerr(stderr);
        clearerr(stdout);

        CHECK(r[0] == FAIL && r[1] == FAIL && r[2] == FAIL && r[3] == OK && r[4] == OK,
              "return values %d %d %d %d %d", r[0], r[1], r[2], r[3], r[4]);

        err = read_file(MSG_ERR_FILE, &n_err);
        out = read_file(MSG_OUT_FILE, &n_out);
        CHECK(out != NULL && n_out == 0, "nothing on stdout (%zu bytes)", n_out);
        CHECK(err != NULL, "stderr captured");
        if(err != NULL){
                const char* expect_head =
                        "samstat: hello world\n"
                        "samstat: warning: w 1\n"
                        "samstat: error: bad value 7 in x.sam\n"
                        "samstat: error: odd value 3\n"
                        "samstat: error: String buffer: ";
                size_t nh = strlen(expect_head);
                CHECK(n_err > nh && strncmp(err, expect_head, nh) == 0, "messages:\n%s", err);
                /* exactly the expected lines: the six short ones (two of them from strbuf_limit),
                 * the long warning and the 100% line */
                int lines = 0;
                for(size_t i = 0; i < n_err; i++){
                        lines += (err[i] == '\n');
                }
                CHECK(lines == 8 && err[n_err - 1] == '\n', "%d lines", lines);
                /* the long warning is complete: prefix, 3000 x, newline */
                char* w = strstr(err, "samstat: warning: xxx");
                CHECK(w != NULL && strspn(w + 18, "x") == 3000 && w[18 + 3000] == '\n', "long message kept whole");
                CHECK(strstr(err, "samstat: 100% done\n") != NULL, "%% in a message");
                CHECK(strstr(err, "Function") == NULL && strstr(err, "failed.") == NULL, "RUN prints nothing");
                CHECK(strstr(err, "n % 2") == NULL, "ASSERT does not print its test");
                CHECK(strstr(err, ".c line") == NULL && strstr(err, "LOG") == NULL && strstr(err, "[20") == NULL,
                      "no source location, type column or time stamp");
                CHECK(strstr(err, "null pointer") == NULL, "MFREE(NULL) / gfree(NULL) are silent");
        }
        free(err);
        free(out);
        remove(MSG_ERR_FILE);
        remove(MSG_OUT_FILE);
        return OK;
}

int main(void)
{
        CHECK(test_messages() == OK, "messages and macros");
        CHECK(test_strbuf() == OK, "string buffers");
        CHECK(test_galloc() == OK, "galloc");
        CHECK(test_files() == OK, "files");
        CHECK(test_revcomp() == OK, "reverse complement");
        CHECK(test_seq_buffer() == OK, "sequence buffers");
        printf("unit_tldlite: %d checks, %d failed\n", n_checks, n_failed);
        return n_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
