#include "sam.h"
#include "tldlite.h"
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define HTSGLUE_IMPORT
#include "htsglue.h"

static void free_aln_data_struct(struct aln_data *a);
static int alloc_aln_data(struct aln_data **aln_d);

int add_aln_data(struct tl_seq_buffer *sb)
{
        ASSERT(sb != NULL,"No sb");
        for(int i = 0; i < sb->malloc_num;i++){
                if(!sb->sequences[i]->data){
                        struct aln_data* f = NULL;
                        RUN(alloc_aln_data(&f));
                        sb->sequences[i]->data = f;
                        f = NULL;
                }
        }
        return OK;
ERROR:
        clear_aln_data(sb);
        return FAIL;
}

int clear_aln_data(struct tl_seq_buffer *sb)
{
        ASSERT(sb != NULL,"No sb");
        for(int i = 0; i < sb->malloc_num;i++){
                if(sb->sequences[i]->data){
                        struct aln_data* a = NULL;

                        a = sb->sequences[i]->data;
                        a->md->len = 0;
                        a->n_clip5 = 0;
                        a->n_clip3 = 0;
                        a->error = -1;
                }
        }
        return OK;
ERROR:
        return FAIL;
}


int free_aln_data(struct tl_seq_buffer *sb)
{
        ASSERT(sb != NULL,"No sb");
        for(int i = 0; i < sb->malloc_num;i++){
                if(sb->sequences[i]->data){
                        free_aln_data_struct(sb->sequences[i]->data);
                        sb->sequences[i]->data = NULL;
                }
        }
        return OK;
ERROR:
        return FAIL;
}


/* After the -p cutoff: does a further primary record follow? Secondary and
 * supplementary records read here lie beyond the cutoff and are not counted. */
static int more_primary_follow(struct sam_bam_file *f)
{
        int r;

        while((r = sam_read1(f->in, f->header, f->b)) >= 0){
                if(f->is_fastx){
                        /* same rule as read_batch: a FASTA/FASTQ record without sequence is not a read */
                        if(f->b->core.l_qseq != 0){
                                return 1;
                        }
                }else if(!(f->b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))){
                        return 1;
                }
        }
        /* an unreadable rest cannot be told apart from more records: report partial */
        return r < -1;
}

/* A slot that held more than this many bytes of record data gives the memory back when
 * the batch is reset: with 4096 slots per batch and several batches in flight, slots must
 * not each keep room for the longest read they ever held. */
#define BATCH_KEEP_DATA 32768

int bam_batch_alloc(struct bam_batch **batch, int cap)
{
        struct bam_batch *b = NULL;

        ASSERT(cap > 0, "A batch needs room for at least one record");
        MMALLOC(b, sizeof(struct bam_batch));
        b->recs = NULL;
        b->n = 0;
        b->alloc = 0;
        b->cap = cap;
        *batch = b;
        return OK;
ERROR:
        return FAIL;
}

void bam_batch_reset(struct bam_batch *batch)
{
        for(int i = 0; i < batch->n; i++){
                bam1_t *r = batch->recs[i];
                if(r->m_data > BATCH_KEEP_DATA){
                        free(r->data);
                        r->data = NULL;
                        r->m_data = 0;
                        r->l_data = 0;
                }
        }
        batch->n = 0;
}

void bam_batch_free(struct bam_batch *batch)
{
        if(batch){
                for(int i = 0; i < batch->alloc; i++){
                        bam_destroy1(batch->recs[i]);
                }
                free(batch->recs);
                MFREE(batch);
        }
}

/* The record slot batch->n, allocated on first use. */
static bam1_t *batch_slot(struct bam_batch *batch)
{
        if(batch->n >= batch->alloc){
                int new_alloc = MACRO_MIN(batch->cap, MACRO_MAX(64, batch->alloc * 2));
                bam1_t **tmp = realloc(batch->recs, sizeof(bam1_t *) * (size_t) new_alloc);
                if(!tmp){
                        return NULL;
                }
                batch->recs = tmp;
                for(int i = batch->alloc; i < new_alloc; i++){
                        batch->recs[i] = bam_init1();
                        if(!batch->recs[i]){
                                batch->alloc = i;
                                return NULL;
                        }
                }
                batch->alloc = new_alloc;
        }
        return batch->recs[batch->n];
}

/* Fill <batch> with at most <max_records> records (the batch limits of the file apply too). */
static int read_batch(struct sam_bam_file *f_handle, struct bam_batch *batch, int max_records)
{
        bam_hdr_t *h = f_handle->header;
        uint64_t n_bases = 0;
        int r = 0;

        batch->n = 0;
        if(f_handle->at_end || f_handle->is_empty || h == NULL){
                /* nothing (more) to read: empty input, EOF or -p cutoff reached */
                return OK;
        }
        if((uint64_t) max_records > f_handle->batch_records){
                max_records = (int) f_handle->batch_records;
        }
        if(max_records > batch->cap){
                max_records = batch->cap;
        }
        ASSERT(max_records > 0, "No room for a record");

        while(1){
                bam1_t *b = batch_slot(batch);
                ASSERT(b != NULL, "Out of memory");
                r = sam_read1(f_handle->in, h, b);
                if(r < 0){
                        break;
                }
                if(b->core.l_qseq < 0 || (uint64_t) b->core.l_qseq > f_handle->max_read_len){
                        ERROR_MSG("%s: record %s has %" PRId32 " bases, more than the %" PRIu64 " a read may have: this is not a file of sequencing reads (a genome FASTA?)",
                                  f_handle->filename, bam_get_qname(b), b->core.l_qseq, f_handle->max_read_len);
                }
                if(!f_handle->is_fastx && (b->core.flag & BAM_FSECONDARY)){
                        f_handle->n_secondary++;
                }else if(!f_handle->is_fastx && (b->core.flag & BAM_FSUPPLEMENTARY)){
                        f_handle->n_supplementary++;
                }else if(f_handle->is_fastx && b->core.l_qseq == 0){
                        /* a FASTA/FASTQ record without sequence is not a read */
                        continue;
                }else{
                        batch->n++;
                        n_bases += (uint64_t) b->core.l_qseq;
                        f_handle->n_primary++;

                        /* -p: stop at exactly max_primary primary records, even mid-batch */
                        if(f_handle->n_primary >= f_handle->max_primary){
                                f_handle->at_end = 1;
                                f_handle->partial = (uint8_t) more_primary_follow(f_handle);
                                return OK;
                        }
                        if(batch->n >= max_records || n_bases >= f_handle->batch_bases){
                                return OK;
                        }
                }
        }
        if(r != -1 && hts_get_format(f_handle->in)->format == cram){
                ERROR_MSG("%s: cannot read the CRAM records (reference sequence not found or file corrupt): set REF_PATH or REF_CACHE to a location holding the reference", f_handle->filename);
        }
        ASSERT(r == -1,"Error reading alignments from %s (truncated or corrupt input?)", f_handle->filename);
        f_handle->at_end = 1;
        return OK;
ERROR:
        return FAIL;
}

int read_bam_batch(struct sam_bam_file *f_handle, struct bam_batch *batch)
{
        return read_batch(f_handle, batch, batch->cap);
}

/* Append len bytes of src to b. Unlike tld_append this does not stop at a 0 byte: the
 * caller gives the length (SEQ and QUAL have one, a 0 inside them must not cut them short).
 * The capacity grows like in tld_append; on failure b is unchanged. */
static int append_bytes(tld_strbuf *b, const char *src, int len)
{
        ASSERT(b != NULL && len >= 0, "Bad append");
        if((int64_t) b->alloc_len - (int64_t) b->len <= (int64_t) len){
                const int64_t want = (int64_t) b->len + (int64_t) len + 1;
                int64_t new_size = MACRO_MAX((int64_t) b->alloc_len * 2, want);
                uint8_t *tmp = NULL;

                ASSERT(want <= INT32_MAX, "String too long");
                if(new_size > INT32_MAX){
                        new_size = INT32_MAX;
                }
                galloc(&tmp, (int) new_size);
                memcpy(tmp, b->str, (size_t) b->len);
                gfree(b->str);
                b->str = tmp;
                b->alloc_len = (int32_t) new_size;
        }
        memcpy(b->str + b->len, src, (size_t) len);
        b->len += len;
        b->str[b->len] = 0;
        return OK;
ERROR:
        return FAIL;
}

/* The largest Phred score the SAM specification allows (QUAL characters '!' to '~'). */
#define PHRED_MAX 93

int convert_bam_batch(const struct bam_batch *batch, int is_fastx, struct tl_seq_buffer *sb)
{
        char* buffer = NULL;
        int b_len = 1024;       /* read lengths are at most max_read_len (read_batch), far below INT_MAX */

        ASSERT(sb->num_seq + batch->n <= sb->malloc_num, "The sequence buffer is too small for the batch");
        galloc(&buffer, b_len);

        for(int k = 0; k < batch->n; k++){
                bam1_t *b = batch->recs[k];
                struct tl_seq* s = NULL;
                struct aln_data* a = NULL;
                uint8_t * seq =    bam_get_seq(b);
                uint8_t* qual_ptr = bam_get_qual(b);

                s = sb->sequences[sb->num_seq];
                a = (struct aln_data* )s->data;

                RUN(tld_append(s->name, bam_get_qname(b)));

                s->len = b->core.l_qseq;

                /* buffer needs room for s->len characters plus the terminating 0 */
                ASSERT(s->len >= 0 && s->len < INT32_MAX - 1024, "Record %s is too long", bam_get_qname(b));
                if(s->len >= b_len){
                        b_len = s->len + 1024;
                        gfree(buffer);
                        buffer = NULL;
                        galloc(&buffer,b_len);
                }
                for (int i = 0; i < s->len; ++i){
                        buffer[i] = seq_nt16_str[bam_seqi(seq, i)];
                }
                buffer[s->len] = 0;

                RUN(append_bytes(s->seq, buffer, s->len));

                if(s->len == 0){
                        /* SEQ '*': no bases and no qualities (qual_ptr must not be read) */
                        RUN(append_bytes(s->qual, buffer, 0));
                }else if(qual_ptr[0] == 0xFF){
                        /* SAM spec: a first quality of 0xFF means QUAL '*' (absent) */
                        if(is_fastx){
                                /* FASTA: no quality string at all */
                                RUN(append_bytes(s->qual, buffer, 0));
                        }else{
                                memset(buffer, 0xFF, s->len);
                                buffer[s->len] = 0;
                                RUN(append_bytes(s->qual, buffer, s->len));
                        }
                }else{
                        /* A BAM can hold any byte as a quality. Above 93 the character
                         * would leave '!'..'~' (and 223 + 33 wraps to 0): clamp, as the
                         * collector would for 94..222 anyway. */
                        for (int i = 0; i < s->len; ++i){
                                buffer[i] = (char) (MACRO_MIN((int) qual_ptr[i], PHRED_MAX) + 33);
                        }
                        buffer[s->len] = 0;
                        RUN(append_bytes(s->qual, buffer, s->len));
                }

                /* FASTA/FASTQ records carry no alignment data (a == NULL). The
                 * paired/read1/read2 flags htslib derives from /1 and /2 name
                 * suffixes are ignored: every read is read 1. */
                if(a){
                        a->flag =  b->core.flag;
                        a->num_hits = 0;
                        a->map_q = b->core.qual;

                        /* NM is stored only when present: reset it so that a reused
                         * slot does not inherit the value of an earlier record. */
                        a->error = -1;
                        if(! (BAM_FUNMAP & b->core.flag)){
                                uint8_t * aux =  bam_aux_get(b,"NM");
                                /* only an integer NM is an edit distance: bam_aux2i reads
                                 * NM:Z / NM:A / NM:f as 0, a perfect match. Anything else,
                                 * or a value that does not fit, is treated as no NM. */
                                if(aux){
                                        switch(aux[0]){
                                        case 'c':
                                        case 'C':
                                        case 's':
                                        case 'S':
                                        case 'i':
                                        case 'I':{
                                                const int64_t nm = bam_aux2i(aux);
                                                if(nm >= 0 && nm <= INT32_MAX){
                                                        a->error = (int) nm;
                                                }
                                                break;
                                        }
                                        default:
                                                break;
                                        }
                                }
                                aux =  bam_aux_get(b,"MD");
                                /* only MD:Z is a mismatch string; bam_aux2Z returns NULL
                                 * for MD:i, MD:A, MD:B ...: treated as no MD */
                                if(aux && aux[0] == 'Z'){
                                        RUN(tld_append(a->md,bam_aux2Z(aux)));
                                }
                        }

                        if((int) b->core.n_cigar > a->n_alloc_cigar){
                                const int new_alloc = MACRO_MAX(a->n_alloc_cigar + a->n_alloc_cigar / 2, (int) b->core.n_cigar);
                                if(a->cigar){
                                        gfree(a->cigar);
                                }
                                a->cigar = NULL;
                                /* nothing allocated until galloc succeeds */
                                a->n_alloc_cigar = 0;
                                a->n_cigar = 0;
                                galloc(&a->cigar, new_alloc);
                                a->n_alloc_cigar = new_alloc;
                        }
                        a->n_cigar = b->core.n_cigar;
                        uint32_t* tmp = bam_get_cigar(b);// ((uint32_t*)((b)->data + (b)->core.l_qname))
                        for(int i = 0; i < a->n_cigar;i++){
                                a->cigar[i] = tmp[i];
                        }

                        a->reverse = bam_is_rev(b);
                }
                if(s->len > sb->max_len){
                        sb->max_len = s->len;
                }
                sb->num_seq++;
        }
        gfree(buffer);
        return OK;
ERROR:
        if(buffer){
                gfree(buffer);
        }
        return FAIL;
}

/* A slot whose read was longer than this gives its memory back (see trim_seq_buffer). */
#define SLOT_KEEP_LEN 32768

static int fresh_strbuf(tld_strbuf **b, int size)
{
        tld_strbuf_free(*b);
        *b = NULL;
        RUN(tld_strbuf_alloc(b, size));
        return OK;
ERROR:
        return FAIL;
}

int trim_seq_buffer(struct tl_seq_buffer *sb, int n_used)
{
        ASSERT(sb != NULL, "No sb");
        ASSERT(n_used >= 0 && n_used <= sb->malloc_num, "Bad slot count");
        for(int i = 0; i < n_used; i++){
                struct tl_seq *s = sb->sequences[i];
                struct aln_data *a = (struct aln_data *) s->data;

                /* tld_strbuf keeps its largest size: replace the buffers of a long read */
                if(s->seq->len > SLOT_KEEP_LEN || s->qual->len > SLOT_KEEP_LEN){
                        RUN(fresh_strbuf(&s->seq, 128));
                        RUN(fresh_strbuf(&s->qual, 128));
                }
                if(a){
                        if(a->md->len > SLOT_KEEP_LEN){
                                RUN(fresh_strbuf(&a->md, 16));
                        }
                        if(a->aln_len_alloc > SLOT_KEEP_LEN){
                                gfree(a->read);
                                gfree(a->genome);
                                a->read = NULL;
                                a->genome = NULL;
                                a->aln_len_alloc = 256;
                                a->aln_len = 0;
                                galloc(&a->read, a->aln_len_alloc);
                                galloc(&a->genome, a->aln_len_alloc);
                        }
                        if(a->n_alloc_cigar > SLOT_KEEP_LEN / 4){
                                gfree(a->cigar);
                                a->cigar = NULL;
                                a->n_alloc_cigar = 16;
                                a->n_cigar = 0;
                                galloc(&a->cigar, a->n_alloc_cigar);
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int read_bam_chunk(struct sam_bam_file *f_handle, struct tl_seq_buffer *sb)
{
        if(!f_handle->chunk){
                int cap = (f_handle->batch_records > (uint64_t) INT32_MAX) ? INT32_MAX : (int) f_handle->batch_records;
                RUN(bam_batch_alloc(&f_handle->chunk, MACRO_MAX(cap, 1)));
        }
        bam_batch_reset(f_handle->chunk);
        RUN(read_batch(f_handle, f_handle->chunk, sb->malloc_num - sb->num_seq));
        RUN(convert_bam_batch(f_handle->chunk, f_handle->is_fastx, sb));
        return OK;
ERROR:
        return FAIL;
}

/* An input that is neither SAM/BAM nor FASTA/FASTQ is acceptable only if it holds nothing
 * but white space after decompression (an empty input). Returns 1 if so, 0 if there is
 * other content (*first is its first character, 0 if the input is more than BLANK_SCAN_MAX
 * bytes of white space), -1 if it cannot be read to the end (a truncated compressed file
 * looks empty to htslib's format detection). The input is read in fixed chunks, never a
 * line at a time: a large file without a newline (binary data, a gzip bomb) must not be
 * loaded into memory. */
#define BLANK_SCAN_MAX (64 * 1024 * 1024)

static int blank_input(samFile *in, int *first)
{
        char buf[65536];
        uint64_t total = 0;

        *first = 0;
        while(total < BLANK_SCAN_MAX){
                ssize_t n;
                if(in->format.compression == no_compression){
                        n = hread(in->fp.hfile, buf, sizeof(buf));
                }else{
                        n = bgzf_read(in->fp.bgzf, buf, sizeof(buf));
                }
                if(n < 0){
                        return -1;
                }
                if(n == 0){
                        return 1;
                }
                for(ssize_t i = 0; i < n; i++){
                        if(!isspace((unsigned char) buf[i])){
                                *first = (unsigned char) buf[i];
                                return 0;
                        }
                }
                total += (uint64_t) n;
        }
        return 0;
}

static void destroy_reader(struct sam_bam_file *f_handle)
{
        if(f_handle->b){
                bam_destroy1(f_handle->b);
        }
        bam_batch_free(f_handle->chunk);
        if(f_handle->header){
                bam_hdr_destroy(f_handle->header);
        }
        if(f_handle->in){
                sam_close(f_handle->in);
        }
        if(f_handle->pool){
                hts_tpool_destroy(f_handle->pool);
        }
        MFREE(f_handle);
}

int open_bam(struct sam_bam_file **file, char *filename)
{
        struct reader_opts o;

        o.pool_threads = 0;
        o.batch_records = 4096;
        o.batch_bases = 8000000;
        o.max_primary = UINT64_MAX;
        o.max_read_len = 0;
        return open_reads(file, filename, &o);
}

/* Opens SAM, BAM, FASTA or FASTQ (plain or compressed; "-" is standard input) and
 * decides what it is from the content, never from the file name. */
int open_reads(struct sam_bam_file **file, char *filename, const struct reader_opts *opts)
{
        struct sam_bam_file* f_handle = NULL;
        const htsFormat *fmt = NULL;
        int r;
        int first = 0;

        MMALLOC(f_handle, sizeof(struct sam_bam_file));

        f_handle->filename = filename;
        f_handle->b = NULL;
        f_handle->chunk = NULL;
        f_handle->in = NULL;
        f_handle->header = NULL;
        f_handle->pool = NULL;
        f_handle->tp.pool = NULL;
        f_handle->tp.qsize = 0;

        f_handle->n_secondary = 0;
        f_handle->n_supplementary = 0;
        f_handle->n_primary = 0;
        f_handle->batch_records = opts->batch_records;
        f_handle->batch_bases = opts->batch_bases;
        f_handle->max_primary = opts->max_primary;
        f_handle->max_read_len = opts->max_read_len ? opts->max_read_len : READER_MAX_READ_LEN;
        f_handle->is_fastx = 0;
        f_handle->is_empty = 0;
        f_handle->has_sq = 0;
        f_handle->at_end = 0;
        f_handle->partial = 0;

        f_handle->in = hts_open(f_handle->filename, "r");
        ASSERT(f_handle->in != NULL, "Cannot open %s (file missing, unreadable or of unknown format)", f_handle->filename);

        fmt = hts_get_format(f_handle->in);
        switch(fmt->format){
        case sam:
        case bam:
                break;
        case cram:
                /* the reference comes from REF_PATH / REF_CACHE only (libcurl is off, so no download) */
                ASSERT(hts_set_opt(f_handle->in, CRAM_OPT_DECODE_MD, 1) == 0, "Cannot enable MD decoding for %s", f_handle->filename);
                break;
        case fasta_format:
        case fastq_format:
                f_handle->is_fastx = 1;
                break;
        case empty_format:
        case text_format:
                /* only white space counts as an (empty) input */
                r = blank_input(f_handle->in, &first);
                ASSERT(r >= 0, "%s: cannot be read (truncated or corrupt compressed file?)", f_handle->filename);
                if(r == 0 && first == '>'){
                        /* htslib takes a FASTA file for nucleotides only if its first sequence
                         * line is made of IUPAC nucleotide letters */
                        ERROR_MSG("%s: looks like FASTA but is not a nucleotide FASTA (protein FASTA is not supported)", f_handle->filename);
                }
                ASSERT(r == 1, "%s: unrecognised input (not SAM, BAM, FASTA or FASTQ)", f_handle->filename);
                f_handle->is_empty = 1;
                break;
        default:
                ERROR_MSG("%s: unrecognised input (not SAM, BAM, FASTA or FASTQ)", f_handle->filename);
                break;
        }

        if(!f_handle->is_empty){
                /* One pool, attached once (htslib >= 1.23.1 refuses a second attach).
                 * BAM: the header is read first. htslib 1.24 deadlocks in bgzf_check_EOF()
                 * when the header of a truncated BAM is read with a pool attached. SAM
                 * and FASTA/FASTQ: before, as a SAM that has read its header takes none. */
                int attach_after_header = (fmt->format == bam);

                if(opts->pool_threads > 0){
                        f_handle->pool = hts_tpool_init(opts->pool_threads);
                        if(!f_handle->pool){
                                WARNING_MSG("Cannot create the decompression thread pool, reading %s without.", f_handle->filename);
                        }
                }
                f_handle->tp.pool = f_handle->pool;
                f_handle->tp.qsize = 0;
                if(f_handle->pool && !attach_after_header){
                        ASSERT(hts_set_opt(f_handle->in, HTS_OPT_THREAD_POOL, &f_handle->tp) == 0, "Cannot attach the thread pool to %s", f_handle->filename);
                }
                f_handle->header = sam_hdr_read(f_handle->in);
                ASSERT(f_handle->header != NULL, "Cannot read the header of %s (truncated or not a SAM/BAM/FASTA/FASTQ file?)", f_handle->filename);
                f_handle->has_sq = sam_hdr_nref(f_handle->header) > 0;
                if(f_handle->pool && attach_after_header){
                        ASSERT(hts_set_opt(f_handle->in, HTS_OPT_THREAD_POOL, &f_handle->tp) == 0, "Cannot attach the thread pool to %s", f_handle->filename);
                }
        }

        f_handle->b = bam_init1();
        ASSERT(f_handle->b != NULL, "bam_init1 failed");

        *file = f_handle;
        return OK;
ERROR:
        if(f_handle){
                destroy_reader(f_handle);
        }
        return FAIL;
}

int close_bam(struct sam_bam_file* f_handle)
{
        int status;

        bam_destroy1(f_handle->b);
        bam_batch_free(f_handle->chunk);
        bam_hdr_destroy(f_handle->header);


        status = sam_close(f_handle->in);
        if(f_handle->pool){
                /* after the file: the pool must outlive every file attached to it */
                hts_tpool_destroy(f_handle->pool);
        }
        if (status < 0) {
                /* report after freeing: the struct must not leak on this path */
                WARNING_MSG("Error closing input: %s.",f_handle->filename);
        }
        MFREE(f_handle);
        ASSERT(status >= 0, "Error closing input");
        return OK;
ERROR:
        return FAIL;
}


int alloc_aln_data(struct aln_data **aln_d)
{
        struct aln_data* a = NULL;
        MMALLOC(a, sizeof(struct aln_data));
        a->md = NULL;
        a->cigar = NULL;
        a->genome = NULL;
        a->read = NULL;
        /* a->clip = NULL; */
        a->aln_len_alloc = 256;
        a->aln_len = 0;
        a->map_q = 0;
        a->num_hits = 0;
        a->n_cigar = 0;
        a->n_alloc_cigar = 16;
        a->reverse = 0;
        a->n_clip5 = 0;
        a->n_clip3 = 0;
        a->error = -1;
        a->m_bases = 0;
        a->i_ev = 0;
        a->i_bases = 0;
        a->d_ev = 0;
        a->d_bases = 0;
        a->s_bases = 0;
        a->x_ops = 0;
        a->has_md_or_x = 0;
        RUN(tld_strbuf_alloc(&a->md, 16));
        galloc(&a->read, a->aln_len_alloc);
        galloc(&a->genome, a->aln_len_alloc);
        galloc(&a->cigar, a->n_alloc_cigar);
        /* galloc(&a->clip, a->aln_len_alloc); */
        *aln_d = a;
        return OK;
ERROR:
        free_aln_data_struct(a);
        return FAIL;
}


void free_aln_data_struct(struct aln_data *a)
{
        if(a){
                if(a->md){
                        tld_strbuf_free(a->md);
                }

                if(a->genome){
                        gfree(a->genome);
                }
                if(a->read){
                        gfree(a->read);
                }

                if(a->cigar){
                        gfree(a->cigar);
                }
                /* if(a->clip){ */
                /*         gfree(a->clip); */
                /* } */
                MFREE(a);
        }
}
