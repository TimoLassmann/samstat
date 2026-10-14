#define REPORT_HELP_IMPORT
#include "report_help.h"

/* Wording follows docs/plan-v2.4/standard-plots.md sections 1.3, 1.4, 4, 5.
   Captions describe shape only; no numeric "typical range" is given.
   Note: bwa-mem and bwa-mem2 write MD tags by default; STAR and minimap2
   (without --MD) do not. Do not claim otherwise in any text here. */

const struct report_plot_help report_plot_help[REPORT_N_PLOTS] = {
        [REPORT_PLOT_LENGTH] = {
                .title = "How long are the reads?",
                .caption = "Reads per length bin (log scale for long reads). The dashed line is N50: half of all bases are in reads at least this long.",
                .good = "One peak where the library preparation or sequencing kit puts it.",
                .look = "A second peak of short reads (adapter dimers, degraded input, failed reads), especially if those reads are unmapped.",
                .benign = "Trimmed short reads form a ramp below the full read length; amplicon libraries show one peak per amplicon.",
                .na_no_sequence = "Not available: the input has no sequence."
        },
        [REPORT_PLOT_READQ] = {
                .title = "How good is each read?",
                .caption = "Reads by their average base quality (higher is better; Q20 = 1 error in 100 bases, Q30 = 1 in 1000).",
                .good = "One peak at high quality.",
                .look = "A second, low-quality peak, especially if it matches the poorly mapped or unmapped reads.",
                .benign = "Instruments that bin qualities give a few spikes instead of a smooth curve.",
                .note = "Base qualities are averaged as error probabilities, so a few bad bases pull a read down and the values read lower than a plain average of Q scores, which is what FastQC reports. Mean Q at the top of the page and plot 3 are averaged in the same way.",
                .na_no_quality = "Not available: the input has no base qualities."
        },
        [REPORT_PLOT_POSQ] = {
                .title = "Where along the read does quality drop?",
                .caption = "Mean base quality at each position along the read, one line per MAPQ class (see the table above).",
                .caption_unaligned = "Mean base quality at each position along the read.",
                .good = "A flat line or a slow decline towards the end of the read.",
                .look = "A sudden drop at one position (a failed cycle) or a steep decline that suggests trimming.",
                .benign = "The first few positions are slightly lower; long-read files use % of read length as the axis, which dilutes effects in the first few bases.",
                .na_no_quality = "Not available: the input has no base qualities."
        },
        [REPORT_PLOT_COMP] = {
                .title = "Is base composition even along the read?",
                .caption = "Percentage of A, C, G, T and N at each position along the read.",
                .good = "Roughly parallel, flat lines.",
                .look = "An N spike at one position, or a change in composition towards the end of the read (in unaligned input often adapter read-through; aligners usually trim adapters off the alignment, which is called soft-clipping, so they are then not shown).",
                .benign = "Uneven first positions in random-primed libraries; amplicon and bisulfite libraries are uneven by design.",
                .na_no_sequence = "Not available: the input has no sequence."
        },
        [REPORT_PLOT_ERRORS] = {
                .title = "Where along the read do errors occur?",
                .caption = "Mismatches, insertions and deletions per aligned base at each position, one line per MAPQ class. Each panel has its own vertical scale.",
                .good = "Low, flat lines, with MAPQ >= 30 the lowest.",
                .look = "Rates that rise towards one end (trim there); errors piled up at the ends of poorly mapped reads (adapters, chimeras); high mismatch rates in well-mapped reads (wrong reference or contamination).",
                .benign = "Low-MAPQ reads are often ambiguous (repeats, reads that fit several places), not necessarily wrong. Long reads carry more insertions and deletions than short reads. A few very long reads can dominate a line.",
                .na_no_mapped = "No mapped reads.",
                .na_no_sequence = "Not available: the input has no sequence."
        },
        [REPORT_PLOT_IDENTITY] = {
                .title = "How accurate are the aligned reads?",
                .caption = "Reads by how closely they match the reference (100% = perfect). The dashed line is the median.",
                .good = "A narrow peak close to 100%.",
                .look = "A broad distribution or a second peak at low identity (contamination, wrong reference, chimeric reads); a high % bases clipped at the top of the page with high identity here (chimeras or adapters hidden by clipping).",
                .benign = "Real differences between the sample and the reference lower identity slightly; short reads pile up at 100%.",
                .note = "Identity = matching bases / aligned columns, each gap counted once; clipped bases are not counted.",
                .na_no_mapped = "No mapped reads.",
                .na_no_data = "Identity not available: no NM tag and no =/X operations."
        }
};

const struct report_mapq_help report_mapq_help[REPORT_N_MAPQ] = {
        [REPORT_MAPQ_GE30]     = { "mapq_ge30",  "MAPQ >= 30", "#0072B2", ""        },
        [REPORT_MAPQ_10_29]    = { "mapq_10_29", "MAPQ 10-29", "#B86E00", "8 4"      },
        [REPORT_MAPQ_LT10]     = { "mapq_lt10",  "MAPQ < 10",  "#D55E00", "2 4"      },
        [REPORT_MAPQ_UNMAPPED] = { "unmapped",   "Unmapped",   "#777777", "12 4 2 4" }
};

const char* const report_how_to_read_label = "How to read plot";   /* + " N" */
const char* const report_note_label        = "How it is calculated:";
const char* const report_good_label        = "Good:";
const char* const report_look_label        = "Worth a look:";
const char* const report_benign_label      = "Often benign:";
const char* const report_mapq_title        = "How well did reads map?";
const char* const report_mapq_caption      = "MAPQ is the aligner's confidence that a read is placed correctly; aligners scale it differently.";

const char* report_plot_caption(int plot_id, int is_aligned)
{
        const struct report_plot_help *h = &report_plot_help[plot_id - 1];
        return (!is_aligned && h->caption_unaligned) ? h->caption_unaligned : h->caption;
}

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
const char* const report_note_no_sequence =
        "Reads have no stored sequence (SEQ is *): length, quality and base composition are not available.";
const char* const report_note_few_reads =
        "Few reads (%s): per-position plots are noisy, and positions covered by fewer than 10 reads are not drawn.";
const char* const report_banner_no_reads =
        "This file contains no reads.";
const char* const report_tile_n50_sub    = "half of all bases are in reads at least this long";
const char* const report_tile_clipped_sub = "trimmed by the aligner (soft-clipped)";
const char* const report_header_partial =
        " \xC2\xB7 partial: first %s reads only (-p); not a random sample";
