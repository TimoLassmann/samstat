#ifndef TLDLITE_H
#define TLDLITE_H

/*
 * tldlite: the few parts of the tld library (github.com/timolassmann/tld) that samstat uses,
 * as one header. It replaces the src/tld git submodule.
 *
 *   - error handling: OK / FAIL, RUN, ASSERT, ERROR_MSG, WARNING_MSG, LOG_MSG
 *   - allocation: MMALLOC, MMALLOC_ALIGNED, MREALLOC, MFREE, galloc / gfree (1D and 2D arrays)
 *   - growable strings: tld_strbuf, tld_append, tld_append_char, TLD_STR, TLD_STRLEN
 *   - files: tld_file_exists, tld_get_filename
 *   - sequence buffers: struct tl_seq, struct tl_seq_buffer, reverse complement
 *
 * Names and behaviour are those of tld, so code written against tld compiles unchanged.
 * tld and samstat have the same author; this copy is distributed under samstat's MIT license
 * (COPYING).
 * Single-header style: every file includes this header; exactly one .c file (tldlite.c)
 * defines TLDLITE_IMPLEMENTATION before including it, which compiles the functions.
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

#define ASSERT(TEST,...)  if(!(TEST)) {                         \
                tld_message(MESSAGE_TYPE_ERROR, #TEST );        \
                tld_message(MESSAGE_TYPE_ERROR, __VA_ARGS__ );  \
                goto ERROR;                                     \
        }

#define LOG_MSG(...) (tld_message(MESSAGE_TYPE_MSG, __VA_ARGS__ ))
#define WARNING_MSG(...) (tld_message(MESSAGE_TYPE_WARNING, __VA_ARGS__ ))
#define ERROR_MSG(...) do {                                     \
                tld_message(MESSAGE_TYPE_ERROR, __VA_ARGS__ );  \
                goto ERROR;                                     \
        }while (0)

#define ADDFAILED(x)  "Function \"" tld_tostring(x) "\" failed."
#define RUN(EXP)                                \
        if((EXP) != OK){                        \
                tld_message(MESSAGE_TYPE_ERROR,ADDFAILED(EXP)); \
                goto ERROR;                                     \
        }

/* One line on stdout: time stamp, LOG / WARNING / ERROR, the message. */
extern void tld_message(uint8_t type, const char *format, ...);

/* ------------------------------------------------------------------ */
/* allocation                                                           */
/* ------------------------------------------------------------------ */

#define MFREE(p)                                        \
        if (p) {                                        \
                free(p);                                \
                (p) = NULL;                             \
        } else {                                        \
                WARNING_MSG("free on a null pointer");  \
        }

#define MMALLOC(p, size)                                        \
        if ((p) != NULL) {                                      \
                ERROR_MSG("malloc on a nun-null pointer");      \
                goto ERROR;                                     \
        }                                                       \
        if ((size) == 0) {                                      \
                ERROR_MSG("malloc of size %d failed", (int)(size)); \
                goto ERROR;                                     \
        }                                                       \
        if (((p) = malloc(size)) == NULL) {                     \
                ERROR_MSG("malloc of size %d failed", (int)(size)); \
                goto ERROR;                                     \
        }                                                       \
        memset(p, 0, size);

/* For a struct declared with more alignment than malloc guarantees (the sequence structs below):
 * plain malloc would break the promise the declaration makes to the compiler. The size is rounded
 * up to a multiple of the alignment (C11 aligned_alloc). Freed with free(), so MFREE works. */
#define MMALLOC_ALIGNED(p, size, align)                                         \
        if ((p) != NULL) {                                                      \
                ERROR_MSG("malloc on a nun-null pointer");                      \
                goto ERROR;                                                     \
        }                                                                       \
        if ((size) == 0) {                                                      \
                ERROR_MSG("malloc of size %d failed", (int)(size));             \
                goto ERROR;                                                     \
        }                                                                       \
        if (((p) = aligned_alloc((align),                                       \
                                 (((size) + (align) - 1) / (align)) * (align))) == NULL) { \
                ERROR_MSG("aligned malloc of size %d (alignment %d) failed",    \
                          (int)(size), (int)(align));                           \
                goto ERROR;                                                     \
        }                                                                       \
        memset(p, 0, size);

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
 * (p[i][j], rows contiguous in one block). gfree(p) frees either. A small header in front of
 * the block records the dimensions. Both jump to ERROR on failure.
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
                *(a) = galloc_hlp(sizeof(__typeof__(**(a))), (b), 0);   \
                if (*(a) == NULL) {                                     \
                        ERROR_MSG("Could not alloc 1D");                \
                }                                                       \
        } while (0);

#define XGALLOC_ARGS_3(a, b, c)                                         \
        do {                                                            \
                if (!tld_valid_2D((a))) {                               \
                        ERROR_MSG("Wrong pointer type for 2D alloc");   \
                }                                                       \
                __typeof__(*(a)) ptr_t = NULL;                          \
                ptr_t = galloc_hlp(sizeof(__typeof__(**(a))), (b), (c)); \
                if (ptr_t == NULL) {                                    \
                        ERROR_MSG("Could not alloc 2D");                \
                }                                                       \
                void *__tld_internal_p = NULL;                          \
                get_ptr(ptr_t, &__tld_internal_p);                      \
                ptr_t[0] = (__typeof__(**(a)))__tld_internal_p;         \
                for (int i = 1; i < (b); i++) {                         \
                        ptr_t[i] = ptr_t[0] + (c)*i;                    \
                }                                                       \
                *(a) = ptr_t;                                           \
        } while (0);

#define XXGALLOC_ARGS(N, ...) XGALLOC_ARGS_##N(__VA_ARGS__)
#define XGALLOC_ARGS(N, ...) XXGALLOC_ARGS(N, __VA_ARGS__)
#define GALLOC_ARGS(...) XGALLOC_ARGS(NARG3(__VA_ARGS__), __VA_ARGS__)
#define galloc(...) GALLOC_ARGS(__VA_ARGS__)
#define gfree(a) tld_free(a)

extern int galloc_too_few(void);
extern void *galloc_hlp(size_t tsize, int dim1, int dim2);
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
} __attribute__((aligned(64)));

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
} __attribute__((aligned(64)));

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
#include <time.h>

/* --- messages --- */

#define TLD_TYPE_MARGIN 8
#define TLD_MESSAGE_MARGIN 22
#define TLD_AT __FILE__ " line " tld_tostring(__LINE__)

static int tld_get_time(char* time_ptr, int size)
{
        struct tm newtime;
        time_t ltime = 0;

        ltime = time(&ltime);
        localtime_r(&ltime, &newtime);
        if(!strftime(time_ptr, size, "[%F %H:%M:%S] ", &newtime)){
                fprintf(stderr,"%*s: ",TLD_TYPE_MARGIN,"ERROR ");
                fprintf(stderr, "Could not get the time");
                fprintf(stderr," (%s)\n",TLD_AT);
                return FAIL;
        }
        return OK;
}

static void tld_vlog(uint8_t type, FILE* f_ptr, const char *format, va_list argp)
{
        char time_string[BUFSIZ];

        if(tld_get_time(time_string,BUFSIZ) != OK){
                fprintf(stderr,"notime");
        }
        fprintf(f_ptr,"%*s: ",TLD_MESSAGE_MARGIN,time_string);
        switch (type) {
        case MESSAGE_TYPE_MSG:
                fprintf(f_ptr,"%*s: ",TLD_TYPE_MARGIN,"LOG ");
                vfprintf(f_ptr, format, argp);
                fprintf(f_ptr,"\n");
                break;
        case MESSAGE_TYPE_WARNING:
                fprintf(f_ptr,"%*s: ",TLD_TYPE_MARGIN,"WARNING ");
                vfprintf(f_ptr, format, argp);
                fprintf(f_ptr,"\n");
                break;
        case MESSAGE_TYPE_ERROR:
                fprintf(f_ptr,"%*s: ",TLD_TYPE_MARGIN,"ERROR ");
                vfprintf(f_ptr, format, argp);
                fprintf(f_ptr," (%s)\n",TLD_AT);
                break;
        default:
                fprintf(f_ptr,"Unrecognized Message type (in tld) - this should never happen");
                break;
        }
        fflush(f_ptr);
}

void tld_message(uint8_t type, const char *format, ...)
{
        va_list argp;
        va_start(argp, format);
        tld_vlog(type,stdout,format, argp);
        va_end(argp);
}

/* --- allocation --- */

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
                tld_mem_i* h = (tld_mem_i*) ((size_t)ptr - tld_mem_shift);
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

void* galloc_hlp(size_t tsize, int dim1, int dim2)
{
        tld_mem_i* h = NULL;
        void* tmp = NULL;

        ASSERT(dim1 >= 1,"DIM1 is too small: %d",dim1);
        MMALLOC(tmp, dim1 * tsize + tld_mem_shift);
        h = (tld_mem_i*)(tmp);
        h->aligned = 0u;
        h->dim1 = dim1;
        h->dim2 = dim2;
        h->ptr = NULL;
        if(dim2){
                MMALLOC(h->ptr, (size_t) dim1 * (size_t) dim2 * tsize);
        }
        return (void *)((size_t)tmp + tld_mem_shift);
ERROR:
        if(tmp){
                free(tmp);
        }
        return NULL;
}

void tld_free(void* p)
{
        tld_mem_i* h = (tld_mem_i*) ((size_t)p - tld_mem_shift);
        if(h->ptr){
                MFREE(h->ptr);
        }
        MFREE(h);
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

static int tld_strbuf_resize(tld_strbuf *b, int new_size)
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

/* room for size more bytes and the terminating 0; capacity at least doubles, so n appends cost O(n) */
static int tld_strbuf_ensure_size(tld_strbuf *b, int size)
{
        int space = b->alloc_len - b->len;
        if(space <= size){
                int new_size = MACRO_MAX(b->alloc_len*2, b->alloc_len + (size - space) + 1);
                RUN(tld_strbuf_resize(b, new_size));
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
        int size = (int) strlen(content);

        RUN(tld_strbuf_ensure_size(b, size));
        memcpy(b->str+b->len,content, size);
        b->len += size;
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
        MMALLOC_ALIGNED(s, sizeof(struct tl_seq), _Alignof(struct tl_seq));
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
        MMALLOC_ALIGNED(sb, sizeof(struct tl_seq_buffer), _Alignof(struct tl_seq_buffer));
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
