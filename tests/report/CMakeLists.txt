# HTML report tests (label: report).
#
# report:write_failure   unwritable output directory gives a clear error
# report:*               content checks on generated HTML (tests/cmake/check_html.cmake)

set(REPORT_HELPERS ${PROJECT_SOURCE_DIR}/tests/cmake)
set(REPORT_OUT ${SAMSTAT_TEST_OUT_DIR}/report)

add_test(NAME report:write_failure
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${SAMSTAT_TEST_DATA_DIR}/xx#minimal.sam"
    -DOUTDIR=${REPORT_OUT}/write_failure
    -P ${REPORT_HELPERS}/run_write_failure.cmake)
set_tests_properties(report:write_failure PROPERTIES LABELS report)

# report_html_test(<name> <input> <check args...>)
# Runs samstat on <input> and applies the given check_html.cmake options
# (-DCONTAINS=... etc.). Option values use '|' as list separator.
function(report_html_test name input)
  add_test(NAME report:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${input}"
      -DOUTDIR=${REPORT_OUT}/${name} ${ARGN}
      -P ${REPORT_HELPERS}/check_html.cmake)
  set_tests_properties(report:${name} PROPERTIES LABELS report)
endfunction()

set(REPORT_SAM "${SAMSTAT_TEST_DATA_DIR}/xx#minimal.sam")
set(REPORT_PAIR_SAM "${SAMSTAT_TEST_DATA_DIR}/xx#pair.sam")
set(REPORT_EMPTY_SAM "${SAMSTAT_FIXTURE_DIR}/empty.sam")

# The page is plain HTML and inline SVG: no script and no Plotly.
foreach(input "xx#minimal.sam" "xx#pair.sam" "ce#1.sam" "ce#2.sam")
  report_html_test(no_script:${input} "${SAMSTAT_TEST_DATA_DIR}/${input}"
    "-DNOT_MATCHES=<script|[Pp]lotly|cdn\\.plot\\.ly"
    "-DCONTAINS=<svg")
endforeach()
report_html_test(no_unused_d3 "${REPORT_SAM}"
  "-DNOT_CONTAINS=d3js.org")

# HTML hygiene
report_html_test(html_lang_viewport "${REPORT_SAM}"
  "-DCONTAINS=<html lang=\"en\">|<meta charset=\"utf-8\">|<meta name=\"viewport\""
  "-DNOT_CONTAINS=100%%")
# A file name with '<' and '&' must be escaped in <title> and <h1> (copy made at test time).
report_html_test(escape_filename "${REPORT_SAM}"
  "-DAS_NAME=a<b&c.sam"
  "-DCONTAINS=<h1>a&lt;b&amp;c.sam</h1>|<title>a&lt;b&amp;c.sam - SAMStat report</title>"
  "-DNOT_CONTAINS=a<b&c.sam</h1>|<title>a<b")

# Footer citations: SAMStat 2 (2023) and the original 2011 paper.
report_html_test(citation "${REPORT_SAM}"
  "-DCONTAINS=10.1093/bioinformatics/btad019|10.1093/bioinformatics/btq614")

# No NaN/inf from divisions by zero (empty input) or anywhere else.
foreach(input "${SAMSTAT_FIXTURE_DIR}/empty.sam" "${SAMSTAT_TEST_DATA_DIR}/xx#minimal.sam" "${SAMSTAT_TEST_DATA_DIR}/xx#pair.sam")
  get_filename_component(iname "${input}" NAME)
  report_html_test(no_nan:${iname} "${input}"
    "-DNOT_MATCHES=NaN|[^A-Za-z]nan[^A-Za-z]|[^A-Za-z]-?inf[^A-Za-z]")
endforeach()

# secondary/supplementary lines are shown only when non-zero
report_html_test(secondary_supplementary_lines "${SAMSTAT_FIXTURE_DIR}/secondary.sam"
  "-DCONTAINS=secondary alignments (not counted)|supplementary alignments (not counted)")
report_html_test(no_secondary_lines_when_zero "${REPORT_SAM}"
  "-DNOT_CONTAINS=secondary alignments (not counted)|supplementary alignments (not counted)")

# =============================================================================
# The report page of 2.4 (WP2.1 wrote the tests, WP4.2a enabled them and added na_notes,
# svg_a11y, no_flags, partial_header and header_line). Specification:
# docs/plan-v2.4/standard-plots.md sections 1, 2 and 5.
# =============================================================================

set(PLOT_T1 "How long are the reads?")
set(PLOT_T2 "How good is each read?")
set(PLOT_T3 "Where along the read does quality drop?")
set(PLOT_T4 "Is base composition even along the read?")
set(PLOT_T5 "Where along the read do errors occur?")
set(PLOT_T6 "How accurate are the aligned reads?")
set(PLOT_MAPQ_T "How well did reads map?")
# aligned input: MAPQ table, then plots 1-6; unaligned input: plots 1-4 and nothing after
set(TITLES_ALIGNED "${PLOT_MAPQ_T}|${PLOT_T1}|${PLOT_T2}|${PLOT_T3}|${PLOT_T4}|${PLOT_T5}|${PLOT_T6}")
set(TITLES_UNALIGNED "${PLOT_T1}|${PLOT_T2}|${PLOT_T3}|${PLOT_T4}")
set(TITLES_ABSENT_UNALIGNED "${PLOT_MAPQ_T}|${PLOT_T5}|${PLOT_T6}")

# plot_set_test(<name> <input> <check args...>): check_plot_set.cmake (titles, no NaN)
function(plot_set_test name input)
  add_test(NAME report:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${input}"
      -DOUTDIR=${REPORT_OUT}/${name} ${ARGN}
      -P ${REPORT_HELPERS}/check_plot_set.cmake)
  set_tests_properties(report:${name} PROPERTIES LABELS report)
endfunction()

# report:plot_set (T5): the titles in order, nothing else; no NaN|Infinity|undefined
plot_set_test(plot_set:short_bam "${SAMSTAT_FIXTURE_DIR}/small.bam" "-DTITLES=${TITLES_ALIGNED}")
plot_set_test(plot_set:ont_like_bam "${SAMSTAT_FIXTURE_DIR}/longread/ont_like.bam" "-DTITLES=${TITLES_ALIGNED}")
plot_set_test(plot_set:long_fq "${SAMSTAT_FIXTURE_DIR}/longread/long.fq"
  "-DTITLES=${TITLES_UNALIGNED}" "-DNOT_TITLES=${TITLES_ABSENT_UNALIGNED}")
# reads of one length: plot 1 is one sentence instead of a histogram of a single bar
plot_set_test(plot_set:single_len_axis "${SAMSTAT_FIXTURE_DIR}/longread/single_len.sam"
  "-DTITLES=${TITLES_ALIGNED}" "-DMATCHES=All [0-9,]+ reads are 150 bp long\\." "-DNOT_CONTAINS=id=\"samstat-p1-t\"|Reads per length bin")

# report:offline: no network access, no external resource (zero requests)
foreach(input small.bam longread/ont_like.bam longread/long.fq)
  string(MAKE_C_IDENTIFIER "${input}" id)
  plot_set_test(offline:${id} "${SAMSTAT_FIXTURE_DIR}/${input}"
    "-DNOT_MATCHES=(src|href)=[\"']?(https?:)?//@@url\\([\"']?https?:")
endforeach()

# report:ultralong_html, report:short_html_size: HTML < 1 MB in both layouts
add_test(NAME report:ultralong_html
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${REPORT_OUT}/ultralong_html -DMODE=ultralong
    -P ${REPORT_HELPERS}/run_html_size.cmake)
add_test(NAME report:short_html_size
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${REPORT_OUT}/short_html_size -DMODE=short
    -P ${REPORT_HELPERS}/run_html_size.cmake)
set_tests_properties(report:ultralong_html report:short_html_size
  PROPERTIES LABELS report)

# =============================================================================
# WP4.2a: what each kind of input shows (standard-plots.md sections 1 and 4)
# =============================================================================

set(LONGREAD "${SAMSTAT_FIXTURE_DIR}/longread")

# Inputs of 49 and 60 reads for the thresholds of the per-position plots (REPORT_MIN_POS_BASES = 50),
# written at configure time into the build tree (no fixture files): single-end 100 bp reads, MAPQ 60,
# qualities that differ per read. mode: perfect = NM:i:0 and MD:Z:100 (error panels all zero),
# mismatch = one mismatch at position i%100 per read, nomd = no NM and no MD tag, twoclass = as
# mismatch, with every odd read at MAPQ 5 (two MAPQ classes of n/2 reads, for the shape view).
set(GEN_DIR ${REPORT_OUT}/gen)
file(MAKE_DIRECTORY ${GEN_DIR})
function(gen_reads path n mode)
  string(REPEAT "ACGT" 25 seq)
  set(txt "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n")
  math(EXPR last "${n}-1")
  foreach(i RANGE 0 ${last})
    math(EXPR pos "1+${i}*10")
    math(EXPR q "33+20+(${i}*7)%20")
    string(ASCII ${q} qc)
    string(REPEAT "${qc}" 100 qual)
    set(tags "")
    set(mapq 60)
    if(mode STREQUAL "twoclass")
      math(EXPR odd "${i}%2")
      if(odd)
        set(mapq 5)
      endif()
    endif()
    if(mode STREQUAL "perfect")
      set(tags "\tNM:i:0\tMD:Z:100")
    elseif(mode STREQUAL "mismatch" OR mode STREQUAL "twoclass")
      math(EXPR mm "${i}%100")
      math(EXPR rest "99-${mm}")
      set(tags "\tNM:i:1\tMD:Z:${mm}A${rest}")
    endif()
    string(APPEND txt "r${i}\t0\tchr1\t${pos}\t${mapq}\t100M\t*\t0\t0\t${seq}\t${qual}${tags}\n")
  endforeach()
  file(WRITE ${path} "${txt}")
endfunction()
gen_reads(${GEN_DIR}/perfect49.sam 49 perfect)
gen_reads(${GEN_DIR}/perfect50.sam 50 perfect)
gen_reads(${GEN_DIR}/perfect60.sam 60 perfect)
gen_reads(${GEN_DIR}/mismatch60.sam 60 mismatch)
gen_reads(${GEN_DIR}/nomd60.sam 60 nomd)
gen_reads(${GEN_DIR}/twoclass120.sam 120 twoclass)
# 60 unaligned reads of 80-139 bp (plots 1-4 are all drawn)
set(fq "")
foreach(i RANGE 0 59)
  math(EXPR len "80+${i}")
  string(REPEAT "A" ${len} s)
  string(REPEAT "I" ${len} q)
  string(APPEND fq "@u${i}\n${s}\n+\n${q}\n")
endforeach()
file(WRITE ${GEN_DIR}/unaligned60.fq "${fq}")
set(NO_PLOT_56 "${PLOT_T5}|${PLOT_T6}|No mapped reads.|Mismatches not available")

# report:shape_view. The key inside each MAPQ-split figure names every class with its read count;
# the collapsed second figure of plots 2 and 6 (plot 1 is the one-length sentence here) scales each
# class to 100 % of its reads. It needs two classes of 50 reads or more, and aligned input.
plot_set_test(shape_view:two_classes "${GEN_DIR}/twoclass120.sam"
  "-DCONTAINS=Compare the MAPQ classes by shape|id=\"samstat-p2s-t\"|id=\"samstat-p6s-t\"|% of the reads in the class|MAPQ &gt;= 30 (60 reads)|MAPQ &lt; 10 (60 reads)"
  "-DNOT_CONTAINS=id=\"samstat-p1s-t\"|id=\"samstat-p3s|id=\"samstat-p5s" "-DSVG_A11Y=1")
plot_set_test(shape_view:one_class "${GEN_DIR}/mismatch60.sam"
  "-DCONTAINS=MAPQ &gt;= 30 (60 reads)" "-DNOT_CONTAINS=Compare the MAPQ classes by shape|% of the reads in the class")
plot_set_test(shape_view:unaligned "${GEN_DIR}/unaligned60.fq"
  "-DNOT_CONTAINS=Compare the MAPQ classes by shape|MAPQ &gt;= 30 (")

# report:na_notes. Unaligned input (FASTQ, FASTA, uBAM without @SQ): plots 1-4 only, no MAPQ
# table, no text of plots 5-6, no MD or identity notes.
plot_set_test(na_notes:fastq "${LONGREAD}/long.fq"
  "-DTITLES=${TITLES_UNALIGNED}" "-DNOT_TITLES=${TITLES_ABSENT_UNALIGNED}"
  "-DNOT_CONTAINS=No mapped reads.|Mismatches not available|Identity|% mapped|MAPQ >= 30"
  "-DSVG_MIN=2")
plot_set_test(na_notes:fasta "${SAMSTAT_FIXTURE_DIR}/reads.fa.gz"
  "-DTITLES=${PLOT_T1}|${PLOT_T4}" "-DNOT_TITLES=${TITLES_ABSENT_UNALIGNED}"
  "-DCONTAINS=Not available: the input has no base qualities."
  "-DNOT_CONTAINS=Mean Q|No mapped reads.")
plot_set_test(na_notes:ubam_like "${LONGREAD}/ubam_like.sam"
  "-DTITLES=${TITLES_UNALIGNED}" "-DNOT_TITLES=${TITLES_ABSENT_UNALIGNED}"
  "-DNOT_CONTAINS=No mapped reads.|Mismatches not available|% mapped|MAPQ >= 30|MAPQ 10-29")
# @SQ present, nothing mapped: aligned layout, the MAPQ table says 100% unmapped, plots 5-6 say so
plot_set_test(na_notes:sq_nomap "${LONGREAD}/sq_nomap.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DMATCHES=<th scope=\"row\">Unmapped</th><td class=\"num\">[0-9,]+</td><td class=\"num\">100\\.0</td>@@No mapped reads\\.</p>")
# X operations stand in for MD: the mismatch row is drawn and no note says reads were left out
plot_set_test(na_notes:eqx_nomd "${LONGREAD}/eqx_nomd.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DNOT_CONTAINS=Mismatches not available|Mismatch profile excludes|Excludes ")
# neither MD nor X: the mismatch row says why, insertions and deletions are still drawn
# (plot 5 keeps its svg), the notes name the reads that were left out (60 reads: plot 5 is drawn)
plot_set_test(na_notes:no_md "${GEN_DIR}/nomd60.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=Mismatches not available: no MD tag.|Mismatch profile excludes 60 of 60 mapped reads|Insertions (% of aligned bases)|Deletions (% of aligned bases)|id=\"samstat-p5-t\"")
# the same with 2 reads: the notes stay, the per-position plots are replaced (see clarity_r2:*)
plot_set_test(na_notes:no_md_few "${SAMSTAT_FIXTURE_DIR}/hardclip_softclip.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=Mismatch profile excludes 2 of 2 mapped reads")
# NM missing: the identity note
plot_set_test(na_notes:nonm "${LONGREAD}/nonm.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=Identity unavailable for ")

# report:svg_a11y: every figure has role="img", <title> and <desc>
# (fixtures of fewer than 50 reads have no per-position plots; mismatch60.sam has all six)
plot_set_test(svg_a11y:aligned "${SAMSTAT_FIXTURE_DIR}/paired_small.sam" "-DSVG_A11Y=1" "-DSVG_MIN=10")
plot_set_test(svg_a11y:ont_like "${LONGREAD}/ont_like.bam" "-DSVG_A11Y=1" "-DSVG_MIN=7")
plot_set_test(svg_a11y:fastq "${LONGREAD}/long.fq" "-DSVG_A11Y=1" "-DSVG_MIN=2")
# (all reads of mismatch60.sam are 100 bp: plot 1 is the sentence, plots 2-6 and the swatches are figures)
plot_set_test(svg_a11y:all_plots "${GEN_DIR}/mismatch60.sam" "-DSVG_A11Y=1" "-DSVG_MIN=10")

# report:no_flags: facts, never verdicts (no PASS, WARN, FAIL; also not as part of a word)
foreach(input small.bam paired_small.sam longread/ont_like.bam longread/long.fq longread/sq_nomap.sam)
  string(MAKE_C_IDENTIFIER "${input}" id)
  plot_set_test(no_flags:${id} "${SAMSTAT_FIXTURE_DIR}/${input}"
    "-DNOT_MATCHES=PASS|WARN|FAIL|[Ww]arning|status")
endforeach()

# report:partial_header: -p puts the partial run in the header line and in the notes
plot_set_test(partial_header "${LONGREAD}/ont_like.sam"
  "-DARGS=-p|5"
  "-DCONTAINS=<span class=\"partial\">|partial: first 5 reads only (-p)|<li>Partial run: first 5 reads only (-p)")
plot_set_test(partial_absent "${LONGREAD}/ont_like.sam"
  "-DNOT_CONTAINS=<span class=\"partial\">|Partial run|(-p)")

# header line (1.1) and summary strip (1.2)
plot_set_test(header_line:aligned_long "${LONGREAD}/ont_like.bam"
  "-DCONTAINS=<h1>ont_like.bam</h1>|BAM &middot|long reads (median 1,434 bp) &middot|SAMStat "
  "-DNOT_CONTAINS=layout|long when"
  "-DTITLES=<dt>Reads</dt>|<dt>Bases</dt>|<dt>N50</dt>|<dt>Mean Q</dt>|<dt>% mapped</dt>|<dt>% MAPQ &gt;= 30</dt>|<dt>% bases clipped</dt>|<dt>Median identity</dt>")
plot_set_test(header_line:aligned_short "${SAMSTAT_FIXTURE_DIR}/paired_small.sam"
  "-DCONTAINS=short reads (median |<dt>Reads (R1+R2)</dt>|<th scope=\"col\" class=\"num\">Read 1</th>"
  "-DNOT_CONTAINS=layout")
plot_set_test(header_line:fastq "${LONGREAD}/long.fq"
  "-DCONTAINS=FASTQ &middot"
  "-DTITLES=<dt>Reads</dt>|<dt>Bases</dt>|<dt>N50</dt>|<dt>Mean Q</dt>"
  "-DNOT_CONTAINS=<dt>% mapped</dt>|<dt>Median identity</dt>")
plot_set_test(header_line:fasta "${SAMSTAT_FIXTURE_DIR}/reads.fa.gz"
  "-DCONTAINS=FASTA, gzip &middot"
  "-DNOT_CONTAINS=<dt>Mean Q</dt>")

# the footer has no command line and no link
plot_set_test(footer:plain "${LONGREAD}/ont_like.bam"
  "-DCONTAINS=<footer>|btad019"
  "-DNOT_CONTAINS=<a href|command line")

# =============================================================================
# Clarity round 1: wording, tiles, captions by input type, exact lengths, degenerate inputs
# (standard-plots.md sections 1.2, 1.3, 2 and 5)
# =============================================================================

# The "What good looks like:" line of every plot is visible under the figure (it is not a verdict on
# this file, so it must not be called "Good:"); the collapsed block names its plot.
plot_set_test(clarity:good_and_howto "${LONGREAD}/ont_like.bam"
  "-DMATCHES=<p class=\"good\"><strong>What good looks like:</strong> One peak where@@<summary>How to read plot 1</summary>@@<summary>How to read plot 6</summary>@@</svg>\n</div>\n<p class=\"good\">"
  "-DNOT_CONTAINS=How to read this|<strong>Good:|<details>\n<summary>How to read plot 1</summary>\n<div class=\"howto\">\n<p><strong>What good looks like:")
# the three percent tiles carry a percent sign; Bases is written in bases
plot_set_test(clarity:tiles "${LONGREAD}/ont_like.bam"
  "-DMATCHES=<dt>% mapped</dt><dd>[0-9.]+%</dd>@@<dt>% MAPQ &gt;= 30</dt><dd>[0-9.]+%</dd>@@<dt>% bases clipped</dt><dd>[0-9.]+%<span class=\"sub\">trimmed by the aligner@@<dt>Bases</dt><dd>89\\.5 kb<span class=\"sub\">89,490</span>@@<dt>N50</dt><dd>[^<]*<span class=\"sub\">half of all bases are in reads at least this long"
  "-DNOT_CONTAINS=89.5 k<")
# the caption of plot 3 names MAPQ classes only for aligned input
# (captions are shown with a figure only: inputs of at least 50 reads, so that plots 3-5 are drawn)
plot_set_test(clarity:caption_aligned "${GEN_DIR}/mismatch60.sam"
  "-DCONTAINS=one line per MAPQ class (see the table above)")
plot_set_test(clarity:caption_unaligned "${GEN_DIR}/unaligned60.fq"
  "-DCONTAINS=Mean base quality at each position along the read.</p>"
  "-DNOT_CONTAINS=MAPQ|mapping-quality")
# plain wording instead of jargon in the visible text; the statistics note is in the collapsed block once
plot_set_test(clarity:wording "${GEN_DIR}/mismatch60.sam"
  "-DCONTAINS=Q20 = 1 error in 100 bases|how closely they match the reference (100% = perfect)|How it is calculated:|Each panel has its own vertical scale."
  "-DNOT_CONTAINS=Phred Q|multimappers|error probabilities, so it reads lower|mapping-quality class|(as in FastQC)")
# jargon that a wet-lab biologist does not know stays out of the page (chimeras, soft-clipping)
plot_set_test(clarity:wording_plain "${GEN_DIR}/mismatch60.sam"
  "-DCONTAINS=A few very long reads can dominate a line.|reads joined from two places|aligners usually cut adapters off the alignment"
  "-DNOT_CONTAINS=himer|oft-clip|soft clip")
# accessibility of the MAPQ table: the class name only in the row header, the table is not named twice
plot_set_test(clarity:mapq_table_a11y "${LONGREAD}/ont_like.bam"
  "-DCONTAINS=<table>|<title id=\"samstat-mapq-0-t\">solid blue line</title>"
  "-DNOT_MATCHES=<table aria-labelledby@@<title id=\"samstat-mapq-[0-9]-t\">MAPQ")
# the page layout: summary values wrap, figures scale to the column, one column on a phone
plot_set_test(clarity:css "${LONGREAD}/ont_like.bam"
  "-DCONTAINS=clamp(1.25rem,5vw,1.75rem)|overflow-wrap:anywhere|@media (max-width:24rem)|.figs svg{display:block;width:100%"
  "-DNOT_CONTAINS=white-space:nowrap")
plot_set_test(clarity:paired_figs "${SAMSTAT_FIXTURE_DIR}/paired_small.sam"
  "-DCONTAINS=<div class=\"figs paired\">")
# few reads: a note (the masking of thin positions is tested in report_plots_test)
plot_set_test(clarity:few_reads "${GEN_DIR}/perfect60.sam"
  "-DCONTAINS=Few reads (60): per-position plots are noisy, and positions covered by fewer than 50 reads are not drawn."
  "-DNOT_CONTAINS=are hidden|Too few reads")
# fewer reads than that: the notice says the plots are hidden, and they are
plot_set_test(clarity:few_reads_hidden "${LONGREAD}/ubam_like.sam"
  "-DCONTAINS=Only 6 reads: per-position plots (3-5) are not reliable and are hidden."
  "-DNOT_CONTAINS=Few reads|fewer than 10")
# an empty file says so and shows nothing else
plot_set_test(clarity:empty_banner "${REPORT_EMPTY_SAM}"
  "-DCONTAINS=<p class=\"banner\">This file contains no reads.</p>|<footer>"
  "-DNOT_CONTAINS=<dl class|n/a|No reads.|<svg|How long are the reads?|How well did reads map?")
# reads without sequence: n/a instead of 0, one note, N/A lines instead of empty axes
plot_set_test(clarity:no_sequence "${SAMSTAT_FIXTURE_DIR}/seq_star_softclip.sam"
  "-DCONTAINS=<dt>Bases</dt><dd>n/a</dd>|<dt>N50</dt><dd>n/a</dd>|Reads have no stored sequence (SEQ is *)|Not available: the input has no sequence."
  "-DNOT_CONTAINS=<dd>0 bp|layout (median|No reads.")

# N50 and median of one 1 Mb read are exact, not a bin edge (run_html_size.cmake MODE=ultralong)
add_test(NAME report:ultralong_exact
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${REPORT_OUT}/ultralong_exact -DMODE=ultralong
    "-DCONTAINS=<dt>N50</dt><dd>1.0 Mb|long reads (median 1.0 Mb)|The read is 1.0 Mb long."
    "-DNOT_CONTAINS=~|1,000.0 kb"
    -P ${REPORT_HELPERS}/run_html_size.cmake)
set_tests_properties(report:ultralong_exact PROPERTIES LABELS report)

# =============================================================================
# Clarity round 2 (reviewers: biologist and accessibility): per-position plots only where they can be
# read, a plot 6 that ends at 100 %, plain subtitle and tile captions, whole-number scales
# =============================================================================

# 49 reads: plots 3-5 are replaced by one line each, the notice and the plots agree, no 10-read rule
plot_set_test(clarity_r2:hidden_below_50 "${GEN_DIR}/perfect49.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=Only 49 reads: per-position plots (3-5) are not reliable and are hidden.|Too few reads (49) for a per-position plot (needs at least 50).|All 49 reads are 100 bp long.|id=\"samstat-p2-t\"|id=\"samstat-p6-t\""
  "-DNOT_CONTAINS=id=\"samstat-p3-t\"|id=\"samstat-p4-t\"|id=\"samstat-p5-t\"|Few reads (|fewer than 10|fewer than 50 reads are not drawn")
# 60 reads: plots 3-5 are drawn; the all-zero error panels say so; the lines carry their colours
plot_set_test(clarity_r2:drawn_from_50 "${GEN_DIR}/perfect60.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=id=\"samstat-p3-t\"|id=\"samstat-p4-t\"|id=\"samstat-p5-t\"|>none detected<|A panel labelled none detected|Positions covered by fewer than 50 reads are not drawn."
  "-DNOT_CONTAINS=Too few reads|are hidden|0.5<"
  "-DSVG_A11Y=1")
# exactly 50 reads: drawn
plot_set_test(clarity_r2:drawn_at_50 "${GEN_DIR}/perfect50.sam"
  "-DCONTAINS=id=\"samstat-p3-t\"|id=\"samstat-p4-t\"|id=\"samstat-p5-t\"|Few reads (50)"
  "-DNOT_CONTAINS=Too few reads|are hidden")
plot_set_test(clarity_r2:errors_drawn "${GEN_DIR}/mismatch60.sam"
  "-DCONTAINS=id=\"samstat-p5-t\""
  "-DNOT_CONTAINS=Too few reads")
# the base colours: contrast >= 5:1 on white, no shared hex with the MAPQ lines (stroke attribute of plot 4)
plot_set_test(clarity_r2:base_colours "${GEN_DIR}/perfect60.sam"
  "-DCONTAINS=stroke=\"#007A5A\"|stroke=\"#A8467F\"|stroke=\"#1F6FA8\""
  "-DNOT_CONTAINS=#CC79A7|#2E86C1|#009E73")
# plot 6 ends at 100: whole-number ticks 95-100 and no tick or bar past 100
plot_set_test(clarity_r2:identity_axis "${GEN_DIR}/perfect60.sam" -DSCOPE=samstat-p6
  "-DCONTAINS=from 95% to 100%, in bins of 1 percentage point|The last bin includes the reads with exactly 100% identity.|>95<|>100<"
  "-DNOT_CONTAINS=>99.2<|>99.4<|>101<")
# subtitle and tile captions in plain words
plot_set_test(clarity_r2:plain_header "${GEN_DIR}/perfect60.sam"
  "-DCONTAINS=short reads (median 100 bp) &middot;|average base quality (Q30 = 1 error in 1000 bases)|reads placed with high confidence|how closely reads match the reference|ends the aligner trimmed off|Files whose median read length is 1 kb or more"
  "-DNOT_CONTAINS=layout|(soft-clipped)")
# one read of 1 Mb: plots 3-5 hidden, plot 1 is a sentence (one length), plot 6 ends at 100
add_test(NAME report:clarity_r2_ultralong
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${REPORT_OUT}/clarity_r2_ultralong -DMODE=ultralong
    "-DCONTAINS=Only 1 read: per-position plots (3-5) are not reliable and are hidden.|Too few reads (1) for a per-position plot (needs at least 50).|The read is 1.0 Mb long.|>95<|>100<"
    "-DNOT_CONTAINS=>1,000.0 kb<|99.2"
    -P ${REPORT_HELPERS}/run_html_size.cmake)
set_tests_properties(report:clarity_r2_ultralong PROPERTIES LABELS report)
