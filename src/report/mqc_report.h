#ifndef MQC_REPORT_H
#define MQC_REPORT_H

#ifdef MQC_REPORT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct stat_collection;
struct samstat_param;

/* Writes MultiQC custom content files next to the HTML report:
     <input>.samstat_general_mqc.json   general statistics columns
     <input>.samstat_mapq_mqc.json      table: reads per MAPQ bin
     <input>.samstat_length_mqc.json    linegraph: read length distribution
     <input>.samstat_basequal_mqc.json  linegraph: mean base quality per position (R1)
   Must be called after stat_collection_finalise(). See docs/multiqc.md. */
EXTERN int mqc_report_write(struct stat_collection* s, struct samstat_param *p, int id);

#undef MQC_REPORT_IMPORT
#undef EXTERN

#endif
