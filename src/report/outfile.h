#ifndef OUTFILE_H
#define OUTFILE_H

#ifdef OUTFILE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#include <stdio.h>

/*
 * Atomic report files, shared by the HTML, JSON and MultiQC writers.
 *
 * outfile_open creates a temporary file .<name>.XXXXXX in the directory of the
 * final name; the caller writes to o->f. outfile_commit flushes, checks for
 * write errors, closes and renames the temporary file over the final name.
 * On any failure (and in outfile_abort) the temporary file is removed and a
 * report that existed before is left as it was. rename() replaces a symlink at
 * the final name instead of following it, so no file outside the report is
 * ever written.
 */
struct outfile {
        FILE *f;
        char *path;             /* final name */
        char *tmp;              /* temporary name in the same directory */
};

/* *path = the formatted string (allocated, any length). */
EXTERN int outfile_path(char **path, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
EXTERN int outfile_open(struct outfile **out, const char *path);
/* Frees *o and sets it to NULL, also on failure. */
EXTERN int outfile_commit(struct outfile **o);
/* Discards the temporary file; NULL-safe, sets *o to NULL. */
EXTERN void outfile_abort(struct outfile **o);

#undef OUTFILE_IMPORT
#undef EXTERN

#endif
