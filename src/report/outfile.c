#include "tldlite.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef HAVE_PTHREADS
#include <pthread.h>
#endif

#define  OUTFILE_IMPORT
#include "outfile.h"

/* mkstemp creates the file with mode 0600; a report gets the mode fopen would
   have given it, 0666 minus the umask. umask() can only be read by setting it,
   and reports are written from several threads, so it is read once. */
static mode_t outfile_umask = 022;

static void outfile_read_umask(void)
{
        outfile_umask = umask(022);
        umask(outfile_umask);
}

static mode_t outfile_mode(void)
{
#ifdef HAVE_PTHREADS
        static pthread_once_t once = PTHREAD_ONCE_INIT;
        pthread_once(&once, outfile_read_umask);
#else
        static int done = 0;
        if(!done){
                outfile_read_umask();
                done = 1;
        }
#endif
        return (mode_t) (0666 & ~outfile_umask);
}

int outfile_path(char **path, const char *fmt, ...)
{
        va_list ap;
        char *buf = NULL;
        int needed;
        int written;

        ASSERT(path != NULL && *path == NULL, "No output for the path");
        va_start(ap, fmt);
        needed = vsnprintf(NULL, 0, fmt, ap);
        va_end(ap);
        ASSERT(needed > 0, "Could not build output file name");
        MMALLOC(buf, (size_t) needed + 1);
        va_start(ap, fmt);
        written = vsnprintf(buf, (size_t) needed + 1, fmt, ap);
        va_end(ap);
        ASSERT(written == needed, "Output file name was truncated");
        *path = buf;
        return OK;
ERROR:
        if(buf){
                free(buf);
        }
        return FAIL;
}

int outfile_open(struct outfile **out, const char *path)
{
        struct outfile *o = NULL;
        struct stat st;
        const char *slash;
        size_t dirlen;
        int fd = -1;

        ASSERT(out != NULL && *out == NULL, "Output file already open");
        MMALLOC(o, sizeof(struct outfile));
        o->f = NULL;
        o->path = NULL;
        o->tmp = NULL;
        o->path = strdup(path);
        ASSERT(o->path != NULL, "out of memory");

        /* A directory at the report name: rename() would fail only after the
           whole report was written; say so now, as fopen did. */
        if(lstat(path, &st) == 0 && S_ISDIR(st.st_mode)){
                ERROR_MSG("Failed to open output file: %s: %s", path, strerror(EISDIR));
        }

        /* the temporary file must be in the same directory: rename() does not
           cross file systems */
        slash = strrchr(path, '/');
        dirlen = slash ? (size_t) (slash - path) + 1 : 0;
        RUN(outfile_path(&o->tmp, "%.*s.%s.XXXXXX", (int) dirlen, path, path + dirlen));

        fd = mkstemp(o->tmp);
        if(fd < 0){
                int e = errno;
                /* nothing was created: there is nothing to unlink */
                free(o->tmp);
                o->tmp = NULL;
                ERROR_MSG("Failed to open output file: %s: %s", path, strerror(e));
        }
        if(fchmod(fd, outfile_mode()) != 0){
                ERROR_MSG("Failed to set the mode of output file: %s: %s", path, strerror(errno));
        }
        o->f = fdopen(fd, "w");
        ASSERT(o->f != NULL, "Failed to open output file: %s: %s", path, strerror(errno));
        fd = -1;                /* now owned by o->f */

        *out = o;
        return OK;
ERROR:
        if(fd >= 0){
                close(fd);
        }
        outfile_abort(&o);
        return FAIL;
}

int outfile_commit(struct outfile **out)
{
        struct outfile *o = NULL;
        int failed;
        int e;

        ASSERT(out != NULL && *out != NULL, "No output file");
        o = *out;
        *out = NULL;

        /* fflush and ferror catch every earlier failed fprintf/fwrite, so the
           writers need not check each call; fclose can still fail (NFS, disk full) */
        errno = 0;
        failed = fflush(o->f) != 0;
        failed |= ferror(o->f) != 0;
        e = errno;
        if(fclose(o->f) != 0 && !failed){
                o->f = NULL;
                ERROR_MSG("Failed to close output file: %s: %s", o->path, strerror(errno));
        }
        o->f = NULL;
        if(failed){
                ERROR_MSG("Failed to write output file: %s: %s", o->path, e ? strerror(e) : "write error");
        }
        if(rename(o->tmp, o->path) != 0){
                ERROR_MSG("Failed to write output file: %s: %s", o->path, strerror(errno));
        }
        free(o->tmp);
        free(o->path);
        free(o);
        return OK;
ERROR:
        outfile_abort(&o);
        return FAIL;
}

void outfile_abort(struct outfile **out)
{
        struct outfile *o;

        if(out == NULL || *out == NULL){
                return;
        }
        o = *out;
        *out = NULL;
        if(o->f){
                fclose(o->f);
        }
        if(o->tmp){
                unlink(o->tmp);
                free(o->tmp);
        }
        if(o->path){
                free(o->path);
        }
        free(o);
}
