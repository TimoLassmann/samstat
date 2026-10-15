#include "tldlite.h"
#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <htslib/bgzf.h>
#include <ctype.h>
#include <inttypes.h>
#include <string.h>
#include <sys/stat.h>

#define PRECHECK_IMPORT
#include "precheck.h"

/* The scanner mirrors how htslib splits the input into records (sam.c, fastq_parse1):
 * FASTA: a record is a '>' line and every following line up to the next '>' line; its
 *        length is the sum of those line lengths.
 * FASTQ: a record is a name line, sequence lines up to a line starting with '+', then
 *        quality lines until they hold as many characters as the sequence (a quality line
 *        may start with '@' or '+').
 * SAM:   lines starting with '@' before the first record are the header; a record is one
 *        line, SEQ is its 10th tab-separated field.
 * Line lengths are counted as htslib does: without the '\n' and without a '\r' before it.
 * The input is read one byte at a time from a fixed buffer; nothing grows with a line. */

/* record names in messages are cut after this many characters */
#define PRECHECK_NAME_MAX 64

enum pc_line_kind {
        LK_NONE = 0,            /* no line started yet */
        LK_HEADER,              /* FASTA '>' line, FASTQ name line */
        LK_SEQ,                 /* FASTA/FASTQ sequence line */
        LK_PLUS,                /* FASTQ '+' line */
        LK_QUAL,                /* FASTQ quality line */
        LK_SAM_HEADER,          /* SAM '@' header line */
        LK_SAM_RECORD,          /* SAM alignment line */
        LK_OTHER                /* an empty SAM line */
};

enum pc_fastq_state {
        FQ_HEADER = 0,          /* the next line is a name line */
        FQ_SEQ,                 /* in the sequence lines */
        FQ_QUAL                 /* in the quality lines */
};

struct pc_scan {
        const char *filename;
        enum htsExactFormat format;
        uint64_t max_read_len;  /* a record may have this many bases */
        uint64_t max_line;      /* a line (other than a SAM record) may have this many characters */
        uint64_t line_no;       /* 1-based number of the current line */
        uint64_t line_len;      /* characters of the current line so far */
        uint64_t seq_len;       /* bases of the current record (SAM: of its SEQ field) so far */
        uint64_t qual_len;      /* FASTQ: quality characters of the current record so far */
        int kind;               /* enum pc_line_kind of the current line */
        int fq_state;           /* enum pc_fastq_state */
        int field;              /* SAM record: 1-based field of the current character */
        int name_len;
        uint8_t at_start;       /* the next character starts a line */
        uint8_t cr;             /* a '\r' was seen and not yet counted */
        uint8_t sam_in_header;  /* no SAM record line seen yet */
        uint8_t name_open;      /* the record name is still being copied */
        uint8_t name_cut;       /* the record name was longer than PRECHECK_NAME_MAX */
        char name[PRECHECK_NAME_MAX + 4];
};

static void name_reset(struct pc_scan *s)
{
        s->name_len = 0;
        s->name[0] = 0;
        s->name_open = 1;
        s->name_cut = 0;
}

static void name_add(struct pc_scan *s, int c)
{
        if(s->name_len >= PRECHECK_NAME_MAX){
                s->name_cut = 1;
                return;
        }
        /* the name goes to the terminal: no control characters */
        s->name[s->name_len++] = (c >= 0x20 && c < 0x7f) ? (char) c : '?';
        s->name[s->name_len] = 0;
}

static const char *name_str(struct pc_scan *s)
{
        if(s->name_len == 0){
                return "(no name)";
        }
        if(s->name_cut){
                /* room for "..." is reserved in the buffer */
                memcpy(s->name + s->name_len, "...", 4);
                s->name_cut = 0;
                s->name_len += 3;
        }
        return s->name;
}

/* SAM: the SEQ field ended (exact) or the scan window ended inside it (not exact). Same
 * wording as the reader's own check (read_batch in htsglue.c), so that the message for a
 * SAM record does not depend on which check finds it first. */
static int sam_seq_done(struct pc_scan *s, int exact)
{
        if(s->seq_len > s->max_read_len){
                ERROR_MSG("%s: record %s has %s%" PRIu64 " bases, more than the %" PRIu64 " a read may have: this is not a file of sequencing reads (a genome FASTA?)",
                          s->filename, name_str(s), exact ? "" : "at least ", s->seq_len, s->max_read_len);
        }
        return OK;
ERROR:
        return FAIL;
}

/* Decide what the line is from its first character (-1: the line is empty). */
static void start_line(struct pc_scan *s, int c)
{
        switch(s->format){
        case fasta_format:
                if(c == '>'){
                        s->kind = LK_HEADER;
                        s->seq_len = 0;
                        name_reset(s);
                }else{
                        s->kind = LK_SEQ;
                }
                break;
        case fastq_format:
                switch(s->fq_state){
                case FQ_HEADER:
                        s->kind = LK_HEADER;
                        s->seq_len = 0;
                        s->qual_len = 0;
                        name_reset(s);
                        s->fq_state = FQ_SEQ;
                        break;
                case FQ_SEQ:
                        if(c == '+'){
                                s->kind = LK_PLUS;
                                s->fq_state = FQ_QUAL;
                        }else{
                                s->kind = LK_SEQ;
                        }
                        break;
                default:
                        s->kind = LK_QUAL;
                        break;
                }
                break;
        default:
                if(c == '@' && s->sam_in_header){
                        s->kind = LK_SAM_HEADER;
                }else if(c < 0){
                        s->kind = LK_OTHER;
                }else{
                        s->sam_in_header = 0;
                        s->kind = LK_SAM_RECORD;
                        s->field = 1;
                        s->seq_len = 0;
                        name_reset(s);
                }
                break;
        }
}

/* One character of the current line ('\n' and a final '\r' are never passed). */
static int add_char(struct pc_scan *s, int c)
{
        if(s->at_start){
                s->at_start = 0;
                start_line(s, c);
        }
        s->line_len++;
        /* a SAM record line holds SEQ, QUAL and tags: only its SEQ is limited (below) */
        if(s->line_len > s->max_line && s->kind != LK_SAM_RECORD){
                ERROR_MSG("%s: line %" PRIu64 " is longer than %" PRIu64 " characters: not a file of sequencing reads (a genome FASTA?)",
                          s->filename, s->line_no, s->max_line);
        }
        switch(s->kind){
        case LK_HEADER:
                /* the name runs from after '>' / '@' to the first white space */
                if(s->line_len > 1 && s->name_open){
                        if(isspace(c)){
                                s->name_open = 0;
                        }else{
                                name_add(s, c);
                        }
                }
                break;
        case LK_SEQ:
                s->seq_len++;
                if(s->seq_len > s->max_read_len){
                        ERROR_MSG("%s: record %s is longer than %" PRIu64 " bases: not a file of sequencing reads (a genome FASTA?)",
                                  s->filename, name_str(s), s->max_read_len);
                }
                break;
        case LK_SAM_RECORD:
                if(c == '\t'){
                        if(s->field == 10){
                                RUN(sam_seq_done(s, 1));
                        }
                        s->field++;
                }else if(s->field == 1){
                        name_add(s, c);
                }else if(s->field == 10){
                        /* checked when the field ends, so that the message gives the
                         * exact length (the field is read to its end, never stored) */
                        s->seq_len++;
                }
                break;
        default:
                break;
        }
        return OK;
ERROR:
        return FAIL;
}

static int end_line(struct pc_scan *s)
{
        if(s->at_start){
                /* an empty line */
                start_line(s, -1);
        }
        if(s->kind == LK_QUAL){
                s->qual_len += s->line_len;
                if(s->qual_len >= s->seq_len){
                        /* htslib reads quality lines until they cover the sequence */
                        s->fq_state = FQ_HEADER;
                }
        }else if(s->kind == LK_SAM_RECORD && s->field == 10){
                RUN(sam_seq_done(s, 1));
        }
        s->line_no++;
        s->line_len = 0;
        s->at_start = 1;
        s->cr = 0;              /* a '\r' before '\n' is not part of the line */
        s->kind = LK_NONE;
        return OK;
ERROR:
        return FAIL;
}

static int scan_chunk(struct pc_scan *s, const char *buf, ssize_t n)
{
        for(ssize_t i = 0; i < n; i++){
                int c = (unsigned char) buf[i];

                if(c == '\n'){
                        RUN(end_line(s));
                        continue;
                }
                if(s->cr){
                        /* the '\r' was not the last character of its line */
                        s->cr = 0;
                        RUN(add_char(s, '\r'));
                }
                if(c == '\r'){
                        s->cr = 1;
                        continue;
                }
                RUN(add_char(s, c));
        }
        return OK;
ERROR:
        return FAIL;
}

int precheck_input(const char *filename, enum htsExactFormat format, uint64_t max_read_len)
{
        char buf[65536];
        struct pc_scan s;
        struct stat st;
        htsFile *in = NULL;
        uint64_t total = 0;

        if(format != sam && format != fasta_format && format != fastq_format){
                return OK;
        }
        /* standard input, a pipe or a process substitution cannot be read twice */
        if(strcmp(filename, "-") == 0 || stat(filename, &st) != 0 || !S_ISREG(st.st_mode)){
                return OK;
        }
        /* a second handle of its own: the reader's handle has only peeked at the content */
        in = hts_open(filename, "r");
        if(in == NULL){
                return OK;
        }
        if(in->format.format != format ||
           !(in->format.compression == no_compression || in->format.compression == gzip || in->format.compression == bgzf)){
                hts_close(in);
                return OK;
        }

        s.filename = filename;
        s.format = format;
        s.max_read_len = max_read_len;
        s.max_line = max_read_len + PRECHECK_LINE_SLACK;
        s.line_no = 1;
        s.line_len = 0;
        s.seq_len = 0;
        s.qual_len = 0;
        s.kind = LK_NONE;
        s.fq_state = FQ_HEADER;
        s.field = 0;
        s.at_start = 1;
        s.cr = 0;
        s.sam_in_header = 1;
        name_reset(&s);
        s.name_open = 0;

        /* only the first PRECHECK_SCAN_MAX bytes: the check must stay cheap for large inputs */
        while(total < PRECHECK_SCAN_MAX){
                size_t want = (size_t) MACRO_MIN((uint64_t) sizeof(buf), (uint64_t) PRECHECK_SCAN_MAX - total);
                ssize_t n;

                if(in->format.compression == no_compression){
                        n = hread(in->fp.hfile, buf, want);
                }else{
                        n = bgzf_read(in->fp.bgzf, buf, want);
                }
                if(n < 0){
                        /* unreadable (truncated or corrupt): left to the reader to report */
                        break;
                }
                if(n == 0){
                        /* end of input: a last line without '\n' ends here */
                        if(!s.at_start || s.cr){
                                s.cr = 0;
                                RUN(end_line(&s));
                        }
                        break;
                }
                RUN(scan_chunk(&s, buf, n));
                total += (uint64_t) n;
        }
        if(total >= PRECHECK_SCAN_MAX && s.kind == LK_SAM_RECORD && s.field == 10){
                /* the scan window ends inside a SEQ field */
                RUN(sam_seq_done(&s, 0));
        }
        hts_close(in);
        return OK;
ERROR:
        if(in){
                hts_close(in);
        }
        return FAIL;
}
