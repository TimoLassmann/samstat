#ifndef SAM_BAM_PARSE_H
#define SAM_BAM_PARSE_H

#include "seq/tld-seq.h"


#ifdef SAM_BAM_PARSE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct  aln_data;

EXTERN int parse_alignment(struct tl_seq* s);

EXTERN int get_readable_cigar(struct aln_data *a, char **cigar);

EXTERN int fix_orientation(struct tl_seq *s);

#undef SAM_BAM_PARSE_IMPORT
#undef EXTERN

#endif
