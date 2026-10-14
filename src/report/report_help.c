#define REPORT_HELP_IMPORT
#include "report_help.h"

/* Wording follows docs/plan-v2.4/standard-plots.md sections 1.3, 1.4, 4, 5.
   Captions describe shape only; no numeric "typical range" is given.
   Note: bwa-mem and bwa-mem2 write MD tags by default; STAR and minimap2
   (without --MD) do not. Do not claim otherwise in any text here. */

const struct report_plot_help report_plot_help[REPORT_N_PLOTS] = {
        [REPORT_PLOT_LENGTH] = {
                "How long are the reads?",
                "Reads per length bin (log scale for long reads). The dashed line is N50: half of all bases are in reads at least this long.",
                "One peak where the library preparation or sequencing kit puts it.",
                "A second peak of short reads (adapter dimers, degraded input, failed reads), especially if those reads are unmapped.",
                "Trimmed short reads form a ramp below the full read length; amplicon libraries show one peak per amplicon.",
                NULL, NULL, NULL
        },
        [REPORT_PLOT_READQ] = {
                "How good is each read?",
                "Reads by their mean base quality. Qualities are averaged as error probabilities, so a few bad bases pull a read down and values read lower than an average of Q scores (as in FastQC).",
                "One peak at high quality.",
                "A second, low-quality peak, especially if it matches the poorly mapped or unmapped reads.",
                "Instruments that bin qualities give a few spikes instead of a smooth curve.",
                "Not available: the input has no base qualities.", NULL, NULL
        },
        [REPORT_PLOT_POSQ] = {
                "Where along the read does quality drop?",
                "Mean base quality at each position along the read, one line per mapping-quality class. Averaged as error probabilities, so it reads lower than FastQC's average of Q scores.",
                "A flat line or a slow decline towards the end of the read.",
                "A sudden drop at one position (a failed cycle) or a steep decline that suggests trimming.",
                "The first few positions are slightly lower; long-read files use % of read length as the axis, which dilutes effects in the first few bases.",
                "Not available: the input has no base qualities.", NULL, NULL
        },
        [REPORT_PLOT_COMP] = {
                "Is base composition even along the read?",
                "Percentage of A, C, G, T and N at each position along the read.",
                "Roughly parallel, flat lines.",
                "An N spike at one position, or a change in composition towards the end of the read (in unaligned input often adapter read-through; aligners usually soft-clip adapters, which are then not shown).",
                "Uneven first positions in random-primed libraries; amplicon and bisulfite libraries are uneven by design.",
                NULL, NULL, NULL
        },
        [REPORT_PLOT_ERRORS] = {
                "Where along the read do errors occur?",
                "Mismatches, insertions and deletions per aligned base at each position, one line per mapping-quality class.",
                "Low, flat lines, with MAPQ >= 30 the lowest.",
                "Rates that rise towards one end (trim there); errors piled up at the ends of poorly mapped reads (adapters, chimeras); high mismatch rates in well-mapped reads (wrong reference or contamination).",
                "Low-MAPQ reads are often ambiguous (repeats, multimappers), not necessarily wrong. Long reads carry more insertions and deletions than short reads.",
                NULL, "No mapped reads.", NULL
        },
        [REPORT_PLOT_IDENTITY] = {
                "How accurate are the aligned reads?",
                "Reads by identity to the reference (matching bases / aligned columns, each gap counted once; clipped bases are not counted). The dashed line is the median.",
                "A narrow peak close to 100%.",
                "A broad distribution or a second peak at low identity (contamination, wrong reference, chimeric reads); a high % bases clipped in the summary with high identity here (chimeras or adapters hidden by clipping).",
                "Real differences between the sample and the reference lower identity slightly; short reads pile up at 100%.",
                NULL, "No mapped reads.",
                "Identity not available: no NM tag and no =/X operations."
        }
};

const struct report_mapq_help report_mapq_help[REPORT_N_MAPQ] = {
        [REPORT_MAPQ_GE30]     = { "mapq_ge30",  "MAPQ >= 30", "#0072B2", ""        },
        [REPORT_MAPQ_10_29]    = { "mapq_10_29", "MAPQ 10-29", "#E69F00", "6 3"     },
        [REPORT_MAPQ_LT10]     = { "mapq_lt10",  "MAPQ < 10",  "#D55E00", "2 3"     },
        [REPORT_MAPQ_UNMAPPED] = { "unmapped",   "Unmapped",   "#777777", "8 3 2 3" }
};

const char* const report_how_to_read_label = "How to read this";
const char* const report_good_label        = "Good:";
const char* const report_look_label        = "Worth a look:";
const char* const report_benign_label      = "Often benign:";
const char* const report_mapq_title        = "How well did reads map?";
const char* const report_mapq_caption      = "MAPQ is the aligner's confidence that a read is placed correctly; aligners scale it differently.";

const char* const report_mismatch_row_na =
        "Mismatches not available: no MD tag. Add one with samtools calmd -b in.bam ref.fa > out.bam (minimap2: --MD).";

const char* const report_note_partial =
        "Partial run: first %s reads only (-p); not a random sample.";
const char* const report_note_no_md =
        "Mismatch profile excludes %s of %s mapped reads (no MD tag and no X operations in CIGAR). Add MD with: samtools calmd -b in.bam ref.fa > out.bam (minimap2: use --MD).";
const char* const report_note_identity_na =
        "Identity unavailable for %s of %s mapped reads (no NM tag and no =/X operations).";
const char* const report_note_beyond_axis =
        "%s%% of bases lie beyond position 1000 and are not in plots 3-5.";
const char* const report_note_mismatch_excludes =
        "Excludes %s of %s mapped reads (no MD tag and no X operations in CIGAR).";
const char* const report_header_partial =
        " \xC2\xB7 partial: first %s reads only (-p); not a random sample";
