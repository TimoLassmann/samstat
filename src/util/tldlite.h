#ifndef TLDLITE_H
#define TLDLITE_H

/*
 * tldlite: the few parts of the tld library (github.com/timolassmann/tld) that samstat uses,
 * as one header. It replaces the src/tld git submodule.
 *
 *   - error handling: OK / FAIL, RUN, ASSERT, ERROR_MSG, WARNING_MSG, LOG_MSG
 *   - allocation: MMALLOC, MREALLOC, MFREE, galloc / gfree (1D and 2D arrays)
 *   - growable strings: tld_strbuf, tld_append, tld_append_char, TLD_STR, TLD_STRLEN
 *   - files: tld_file_exists, tld_get_filename
 *   - sequence buffers: struct tl_seq, struct tl_seq_buffer, reverse complement
 *
 * Names follow tld, so code written against tld compiles unchanged. Unlike tld, messages go to
 * stderr as "samstat: error: ..." without time stamps or source locations, and RUN prints
 * nothing (the function that failed has already said why).
 * tld and samstat have the same author; this copy is distributed under samstat's MIT license
 * (COPYING).
 * Single-header style: every file includes this header; exactly one .c file (tldlite.c)
 * defines TLDLITE_IMPLEMENTATION before including it, which compiles the functions.
 *
 * Every statement-like macro (RUN, ASSERT, ERROR_MSG, MFREE, MMALLOC, MREALLOC, galloc) is a
 * single statement written as do { ... } while (0): it must be followed by ';' and is safe as
 * the body of an if / else without braces.
 */

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* error handling                                                       */
/* ------------------------------------------------------------------ */

#define tld_stringify(x) #x
#define tld_tostring(x) tld_stringify(x)

#define OK              0
#define FAIL            1

#define MACRO_MIN(a,b)          (((a)<(b))?(a):(b))
#define MACRO_MAX(a,b)          (((a)>(b))?(a):(b))

#define MESSAGE_TYPE_MSG 0
#define MESSAGE_TYPE_WARNING 1
#define MESSAGE_TYPE_ERROR 2

/* Only the message is printed: the stringified test is code, not a message, and passing it
 * as a format would be undefined behaviour for a test such as (n % 2 == 0). */
#define ASSERT(TEST,...) do {                                   \
                if(!(TEST)) {                                   \
                        tld_message(MESSAGE_TYPE_ERROR, __VA_ARGS__ ); \
                        goto ERROR;                             \
                }                                               \
        } while (0)

#define LOG_MSG(...) (tld_message(MESSAGE_TYPE_MSG, __VA_ARGS__ ))
#define WARNING_MSG(...) (tld_message(MESSAGE_TYPE_WARNING, __VA_ARGS__ ))
#define ERROR_MSG(...) do {                                     \
                tld_message(MESSAGE_TYPE_ERROR, __VA_ARGS__ );  \
                goto ERROR;                                     \
        } while (0)

/* RUN prints nothing: the function that failed has printed the cause with ERROR_MSG / ASSERT.
 * (tld printed "Function \"<code>\" failed." with the code text as the printf format.) */
#define RUN(EXP) do {                                           \
                if((EXP) != OK){                                \
                        goto ERROR;                             \
                }                                               \
        } while (0)

/* One line on stderr, written in one piece: "samstat: <msg>" (LOG), "samstat: warning: <msg>"
 * or "samstat: error: <msg>". */
extern void tld_message(uint8_t type, const char *format, ...) __attribute__((format(printf, 2, 3)));

/* ------------------------------------------------------------------ */
/* allocation                                                           */
/* ------------------------------------------------------------------ */

/* Frees p and sets it to NULL; a NULL p is a no-op, as for free(). */
#define MFREE(p) do {                                           \
                if (p) {                                        \
                        free(p);                                \
                        (p) = NULL;                             \
                }                                               \
        } while (0)

/* p must be NULL; on success it points to size zeroed bytes. size is evaluated once. */
#define MMALLOC(p, size) do {                                   \
                size_t tld_mm_size_ = (size_t) (size);          \
                if ((p) != NULL) {                              \
                        ERROR_MSG("malloc on a non-null pointer"); \
                }                                               \
                if (tld_mm_size_ == 0) {                        \
                        ERROR_MSG("malloc of size 0");          \
                }                                               \
                if (((p) = malloc(tld_mm_size_)) == NULL) {     \
                        ERROR_MSG("malloc of size %zu failed", tld_mm_size_); \
                }                                               \
                memset((p), 0, tld_mm_size_);                   \
        } while (0)

#define MREALLOC(p, size)                                               \
        do {                                                            \
                void *tmpp;                                             \
                if ((size) == 0) {                                      \
                        ERROR_MSG("malloc of size %d failed", (int)(size)); \
                        goto ERROR;                                     \
                }                                                       \
                if ((p) == NULL) {                                      \
                        tmpp = malloc(size);                            \
                } else {                                                \
                        tmpp = realloc((p), (size));                    \
                }                                                       \
                if (tmpp != NULL) {                                     \
                        (p) = tmpp;                                     \
                } else {                                                \
                        ERROR_MSG("realloc for size %d failed", (int)(size)); \
                        goto ERROR;                                     \
                }                                                       \
        } while (0)

/*
 * galloc(&p, n) allocates n zeroed elements of *p's type; galloc(&p, n, m) an n x m array
 * (p[i][j], rows contiguous in one block). gfree(p) frees either; gfree(NULL) is a no-op.
 * A small header in front of the block records the dimensions. Both jump to ERROR on failure.
 */
#define tld_valid_1D(x)                                                 \
        _Generic((x), char ** : 1, int8_t ** : 1, uint8_t ** : 1, int16_t ** : 1, \
                 uint16_t ** : 1, int32_t ** : 1, uint32_t ** : 1, int64_t ** : 1, \
                 uint64_t ** : 1, float ** : 1, double ** : 1, default : 0)

#define tld_valid_2D(x)                                                 \
        _Generic((x), char *** : 1, int8_t *** : 1, uint8_t *** : 1,    \
                 int16_t *** : 1, uint16_t *** : 1, int32_t *** : 1,    \
                 uint32_t *** : 1, int64_t *** : 1, uint64_t *** : 1, float *** : 1, \
                 double *** : 1, default : 0)

#define XARG3(_0, _1, _2, _3, ...) _3
#define NARG3(...) XARG3(__VA_ARGS__, 3, 2, 1, 0)

#define XGALLOC_ARGS_1(a) galloc_too_few()
#define XGALLOC_ARGS_2(a, b)                                            \
        do {                                                            \
                if (!tld_valid_1D((a))) {                               \
                        ERROR_MSG("Wrong pointer type for 1D alloc");   \
                }                                                       \
                *(a) = galloc_hlp(sizeof(**(a)), 0, (b), 0);            \
                if (*(a) == NULL) {                                     \
                        ERROR_MSG("Could not alloc 1D");                \
                }                                                       \
        } while (0)

/* The row pointers have the size sizeof(**a); the data block holds b * c elements of
 * sizeof(***a) (sizing it with the pointer size overflows for 8-byte elements on 32-bit). */
#define XGALLOC_ARGS_3(a, b, c)                                         \
        do {                                                            \
                if (!tld_valid_2D((a))) {                               \
                        ERROR_MSG("Wrong pointer type for 2D alloc");   \
                }                                                       \
                __typeof__(*(a)) tld_galloc_rows_ = NULL;               \
                tld_galloc_rows_ = galloc_hlp(sizeof(**(a)), sizeof(***(a)), (b), (c)); \
                if (tld_galloc_rows_ == NULL) {                         \
                        ERROR_MSG("Could not alloc 2D");                \
                }                                                       \
                void *tld_galloc_data_ = NULL;                          \
                get_ptr(tld_galloc_rows_, &tld_galloc_data_);           \
                tld_galloc_rows_[0] = (__typeof__(**(a))) tld_galloc_data_; \
                for (int tld_galloc_i_ = 1; tld_galloc_i_ < (b); tld_galloc_i_++) { \
                        tld_galloc_rows_[tld_galloc_i_] = tld_galloc_rows_[0] + (size_t) (c) * (size_t) tld_galloc_i_; \
                }                                                       \
                *(a) = tld_galloc_rows_;                                \
        } while (0)

#define XXGALLOC_ARGS(N, ...) XGALLOC_ARGS_##N(__VA_ARGS__)
#define XGALLOC_ARGS(N, ...) XXGALLOC_ARGS(N, __VA_ARGS__)
#define GALLOC_ARGS(...) XGALLOC_ARGS(NARG3(__VA_ARGS__), __VA_ARGS__)
#define galloc(...) GALLOC_ARGS(__VA_ARGS__)
#define gfree(a) tld_free(a)

extern int galloc_too_few(void);
/* dim1 elements of tsize bytes; with dim2 > 0 also a data block of dim1 * dim2 elements of
 * esize bytes (reached with get_ptr). NULL on failure. */
extern void *galloc_hlp(size_t tsize, size_t esize, int dim1, int dim2);
extern void tld_free(void *p);
extern int get_ptr(void *ptr, void **d);

/* ------------------------------------------------------------------ */
/* growable strings                                                     */
/* ------------------------------------------------------------------ */

/* str[len] is always 0; alloc_len is the capacity of str. */
typedef struct tld_string_buffer {
        uint8_t* str;
        int32_t len;
        int32_t alloc_len;
} tld_strbuf;

#define TLD_STR(X) ((char*)(X)->str)
#define TLD_STRLEN(X) (X)->len

extern int tld_strbuf_alloc(tld_strbuf **buffer, int size);   /* size 0: 1024 */
extern void tld_strbuf_free(tld_strbuf *b);
extern int tld_strbuf_clear(tld_strbuf *b);
extern int tld_append(tld_strbuf *b, char *content);
extern int tld_append_char(tld_strbuf *b, char c);

/* ------------------------------------------------------------------ */
/* files                                                                */
/* ------------------------------------------------------------------ */

/* OK if the path exists. */
extern int tld_file_exists(const char* name);
/* *out = a malloc'ed copy of the part of path after the last '/' (free with MFREE). */
extern int tld_get_filename(char* path, char** out);

/* ------------------------------------------------------------------ */
/* sequence buffers                                                     */
/* ------------------------------------------------------------------ */

#define TL_SEQ_BUFFER_DNA 4

struct tl_seq{
        tld_strbuf* seq;
        tld_strbuf* name;
        tld_strbuf* qual;
        void* data;
        int len;
};

struct tl_seq_buffer{
        struct tl_seq** sequences;
        void* data;
        int malloc_num;
        int offset;
        int base_quality_offset;
        int num_seq;
        int max_len;
        int is_fastq;
        int L;
};

extern int alloc_tl_seq_buffer(struct tl_seq_buffer** seq_buf, int size);
extern int reset_tl_seq_buffer(struct tl_seq_buffer* sb);      /* every length to 0, num_seq to 0 */
extern void free_tl_seq_buffer(struct tl_seq_buffer* sb);
extern int alloc_tl_seq(struct tl_seq** sequence);
extern void free_tl_seq(struct tl_seq* sequence);

/* In place: reverse l bytes; reverse-complement l bases (ACGT and acgt; other bytes unchanged). */
extern int reverse(uint8_t *s, int l);
extern int reverse_comp(uint8_t *s, int l);
/* seq reverse-complemented, qual reversed */
extern int rev_comp_tl_seq(struct tl_seq* s);

#endif /* TLDLITE_H */

/* ================================================================== */
/* implementation                                                       */
/* ================================================================== */

#ifdef TLDLITE_IMPLEMENTATION
#ifndef TLDLITE_IMPLEMENTATION_DONE
#define TLDLITE_IMPLEMENTATION_DONE

#include <sys/stat.h>

/* --- messages --- */

/* Messages from worker threads must not interleave: the whole line is formatted first and
 * written with one fwrite under the stream lock. Most lines fit the stack buffer; a longer
 * one is formatted again into a heap buffer (and truncated if even that fails). */
#define TLD_MESSAGE_BUF 1024

static void tld_vlog(uint8_t type, FILE* f_ptr, const char *format, va_list argp)
{
        char stack_buf[TLD_MESSAGE_BUF] = {0};
        char* heap_buf = NULL;
        char* line = stack_buf;
        const char* prefix = NULL;
        size_t prefix_len = 0;
        size_t len = 0;
        int n = 0;
        va_list again;

        switch (type) {
        case MESSAGE_TYPE_WARNING:
                prefix = "samstat: warning: ";
                break;
        case MESSAGE_TYPE_ERROR:
                prefix = "samstat: error: ";
                break;
        default:
                prefix = "samstat: ";
                break;
        }
        prefix_len = strlen(prefix);
        memcpy(stack_buf, prefix, prefix_len);

        va_copy(again, argp);
        /* room for the newline and the terminating 0 */
        n = vsnprintf(stack_buf + prefix_len, sizeof(stack_buf) - prefix_len - 1, format, argp);
        if(n < 0){
                static const char bad[] = "(message could not be formatted)";
                memcpy(stack_buf + prefix_len, bad, sizeof(bad));
                len = prefix_len + sizeof(bad) - 1;
        }else if((size_t) n < sizeof(stack_buf) - prefix_len - 1){
                len = prefix_len + (size_t) n;
        }else{
                size_t need = prefix_len + (size_t) n + 2;
                heap_buf = malloc(need);
                if(heap_buf != NULL){
                        memcpy(heap_buf, prefix, prefix_len);
                        vsnprintf(heap_buf + prefix_len, need - prefix_len - 1, format, again);
                        line = heap_buf;
                        len = prefix_len + (size_t) n;
                }else{
                        len = sizeof(stack_buf) - 2;   /* the truncated message */
                }
        }
        va_end(again);
        line[len] = '\n';
        line[len + 1] = 0;

        flockfile(f_ptr);
        fwrite(line, 1, len + 1, f_ptr);
        fflush(f_ptr);
        funlockfile(f_ptr);
        free(heap_buf);
}

void tld_message(uint8_t type, const char *format, ...)
{
        va_list argp;
        va_start(argp, format);
        tld_vlog(type, stderr, format, argp);
        va_end(argp);
}

/* --- allocation --- */

/* aligned(16): the header is 16 bytes on 32- and 64-bit, so the elements after it keep
 * malloc's alignment */
typedef struct {
        uint32_t aligned : 1;
        uint32_t dim1 : 31;
        uint32_t dim2;
        void* ptr;
} __attribute__((aligned(16))) tld_mem_i;

#define tld_mem_shift sizeof(tld_mem_i)

int get_ptr(void* ptr, void** d)
{
        if(ptr){
                tld_mem_i* h = (tld_mem_i*) (void*) ((char*) ptr - tld_mem_shift);
                *d = h->ptr;
                return OK;
        }
        return FAIL;
}

int galloc_too_few(void)
{
        ERROR_MSG("Too few parameters passed to galloc.\nExpects at least 2");
        return OK;
ERROR:
        return FAIL;
}

void* galloc_hlp(size_t tsize, size_t esize, int dim1, int dim2)
{
        tld_mem_i* h = NULL;
        void* tmp = NULL;

        ASSERT(dim1 >= 1,"DIM1 is too small: %d",dim1);
        ASSERT(dim2 >= 0,"DIM2 is negative: %d",dim2);
        ASSERT(tsize > 0 && (size_t) dim1 <= (SIZE_MAX - tld_mem_shift) / tsize,
               "galloc: %d elements of %zu bytes do not fit in memory", dim1, tsize);
        if(dim2){
                ASSERT(esize > 0 && (size_t) dim2 <= SIZE_MAX / esize / (size_t) dim1,
                       "galloc: %d x %d elements of %zu bytes do not fit in memory", dim1, dim2, esize);
        }
        MMALLOC(tmp, (size_t) dim1 * tsize + tld_mem_shift);
        h = (tld_mem_i*)(tmp);
        h->aligned = 0u;
        h->dim1 = (uint32_t) dim1;
        h->dim2 = (uint32_t) dim2;
        h->ptr = NULL;
        if(dim2){
                MMALLOC(h->ptr, (size_t) dim1 * (size_t) dim2 * esize);
        }
        return (void *)((char*) tmp + tld_mem_shift);
ERROR:
        if(tmp){
                free(tmp);
        }
        return NULL;
}

/* gfree(NULL) is a no-op, as for free() */
void tld_free(void* p)
{
        tld_mem_i* h = NULL;
        if(p == NULL){
                return;
        }
        h = (tld_mem_i*) (void*) ((char*) p - tld_mem_shift);
        free(h->ptr);
        free(h);
}

/* --- growable strings --- */

int tld_strbuf_alloc(tld_strbuf **buffer, int size)
{
        tld_strbuf* b = NULL;
        MMALLOC(b,sizeof(tld_strbuf));
        if(size == 0){
                size = 1024;
        }
        b->alloc_len = size;
        b->len = 0;
        b->str = NULL;
        galloc(&b->str, b->alloc_len);
        *buffer = b;
        return OK;
ERROR:
        if(b){
                MFREE(b);
        }
        return FAIL;
}

void tld_strbuf_free(tld_strbuf *b)
{
        if(b){
                if(b->str){
                        gfree(b->str);
                }
                MFREE(b);
        }
}

static int tld_strbuf_resize(tld_strbuf *b, int32_t new_size)
{
        uint8_t* tmp = NULL;
        if(new_size > b->alloc_len){
                galloc(&tmp, new_size);
                memcpy(tmp,b->str,b->len);
                tmp[b->len] = 0;
                gfree(b->str);
                b->str = tmp;
                b->alloc_len = new_size;
        }
        return OK;
ERROR:
        if(tmp){
                gfree(tmp);
        }
        return FAIL;
}

/* Room for size more bytes and the terminating 0; capacity at least doubles, so n appends
 * cost O(n). The arithmetic is in size_t: len and alloc_len are int32_t, so alloc_len * 2
 * overflowed in int above 1 GiB. A buffer is capped at INT32_MAX bytes; more fails cleanly. */
static int tld_strbuf_ensure_size(tld_strbuf *b, size_t size)
{
        size_t len = (size_t) b->len;
        size_t cap = (size_t) b->alloc_len;

        if(cap - len <= size){
                ASSERT(size < (size_t) INT32_MAX - len,
                       "String buffer: %zu more bytes do not fit (%zu used, at most %d)", size, len, INT32_MAX);
                size_t new_size = MACRO_MAX(cap * 2, len + size + 1);
                new_size = MACRO_MIN(new_size, (size_t) INT32_MAX);
                RUN(tld_strbuf_resize(b, (int32_t) new_size));
        }
        return OK;
ERROR:
        return FAIL;
}

int tld_strbuf_clear(tld_strbuf *b)
{
        ASSERT(b != NULL,"No buffer");
        b->len = 0;
        b->str[b->len] = 0;
        return OK;
ERROR:
        return FAIL;
}

int tld_append(tld_strbuf *b, char *content)
{
        size_t size = strlen(content);

        /* fails if b would grow past INT32_MAX bytes, so the int32_t len cannot overflow */
        RUN(tld_strbuf_ensure_size(b, size));
        memcpy(b->str+b->len,content, size);
        b->len += (int32_t) size;
        b->str[b->len] = 0;
        return OK;
ERROR:
        return FAIL;
}

int tld_append_char(tld_strbuf *b, char c)
{
        RUN(tld_strbuf_ensure_size(b, 1));
        b->str[b->len] = c;
        b->len++;
        b->str[b->len] = 0;
        return OK;
ERROR:
        return FAIL;
}

/* --- files --- */

int tld_file_exists(const char* name)
{
        struct stat buf;
        if(stat(name, &buf) == 0){
                return OK;
        }
        return FAIL;
}

int tld_get_filename(char* path, char** out)
{
        char* tmp = NULL;
        int len = (int) strlen(path);
        int c = 0;

        MMALLOC(tmp, sizeof(char) * (len+1));
        for(int i = 0; i < len; i++){
                tmp[c] = path[i];
                c++;
                if(path[i] == '/'){
                        c = 0;
                }
        }
        tmp[c] = 0;
        if(c == 0){
                ERROR_MSG("No filename found in: %s", path);
        }
        *out = tmp;
        return OK;
ERROR:
        if(tmp){
                MFREE(tmp);
        }
        return FAIL;
}

/* --- sequence buffers --- */

int alloc_tl_seq(struct tl_seq** sequence)
{
        struct tl_seq* s = NULL;
        MMALLOC(s, sizeof(struct tl_seq));
        s->len = 0;
        s->seq = NULL;
        s->qual = NULL;
        s->name = NULL;
        s->data = NULL;
        RUN(tld_strbuf_alloc(&s->seq, 128));
        RUN(tld_strbuf_alloc(&s->qual, 128));
        RUN(tld_strbuf_alloc(&s->name, 128));
        *sequence = s;
        return OK;
ERROR:
        free_tl_seq(s);
        return FAIL;
}

void free_tl_seq(struct tl_seq* sequence)
{
        if(sequence){
                if(sequence->seq){
                        tld_strbuf_free(sequence->seq);
                }
                if(sequence->qual){
                        tld_strbuf_free(sequence->qual);
                }
                if(sequence->name){
                        tld_strbuf_free(sequence->name);
                }
                MFREE(sequence);
        }
}

int alloc_tl_seq_buffer(struct tl_seq_buffer** seq_buf, int size)
{
        struct tl_seq_buffer* sb = NULL;
        ASSERT(size > 0, "Size of sequence buffer is %d",size);
        MMALLOC(sb, sizeof(struct tl_seq_buffer));
        sb->malloc_num = size;
        sb->max_len = 0;
        sb->num_seq = 0;
        sb->L = 0;
        sb->is_fastq = 0;
        sb->offset = 0;
        sb->base_quality_offset = 0;
        sb->sequences = NULL;
        sb->data = NULL;
        MMALLOC(sb->sequences, sizeof(struct tl_seq*) * sb->malloc_num);
        for(int i = 0; i < sb->malloc_num; i++){
                sb->sequences[i] = NULL;
                RUN(alloc_tl_seq(&sb->sequences[i]));
        }
        *seq_buf = sb;
        return OK;
ERROR:
        free_tl_seq_buffer(sb);
        return FAIL;
}

int reset_tl_seq_buffer(struct tl_seq_buffer* sb)
{
        ASSERT(sb != NULL, "No sequence buffer");
        for(int i = 0; i < sb->malloc_num; i++){
                sb->sequences[i]->len = 0;
                sb->sequences[i]->name->len = 0;
                sb->sequences[i]->seq->len = 0;
                sb->sequences[i]->qual->len = 0;
        }
        sb->num_seq = 0;
        return OK;
ERROR:
        return FAIL;
}

void free_tl_seq_buffer(struct tl_seq_buffer* sb)
{
        if(sb){
                if(sb->sequences){
                        for(int i = 0; i < sb->malloc_num; i++){
                                free_tl_seq(sb->sequences[i]);
                        }
                        MFREE(sb->sequences);
                }
                MFREE(sb);
        }
}

static uint8_t tld_comp(uint8_t c)
{
        switch(c){
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default: return c;
        }
}

int reverse(uint8_t *s, int l)
{
        for(int first = 0, last = l - 1; first < last; first++, last--){
                uint8_t tmp = s[first];
                s[first] = s[last];
                s[last] = tmp;
        }
        return OK;
}

int reverse_comp(uint8_t *s, int l)
{
        int first = 0;
        int last = l - 1;
        for(; first < last; first++, last--){
                uint8_t tmp = tld_comp(s[first]);
                s[first] = tld_comp(s[last]);
                s[last] = tmp;
        }
        if(first == last){
                s[first] = tld_comp(s[first]);     /* the middle base of an odd length */
        }
        return OK;
}

int rev_comp_tl_seq(struct tl_seq* s)
{
        if(s->seq){
                reverse_comp(s->seq->str, s->seq->len);
        }
        if(s->qual){
                reverse(s->qual->str, s->qual->len);
        }
        return OK;
}

#endif /* TLDLITE_IMPLEMENTATION_DONE */
#endif /* TLDLITE_IMPLEMENTATION */
