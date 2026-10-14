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
# the short layout x axis ends at the longest read (150 bp for single_len.sam), not at 1000.
# (the tick label is the text node >150< of the x axis)
plot_set_test(plot_set:single_len_axis "${SAMSTAT_FIXTURE_DIR}/longread/single_len.sam"
  "-DTITLES=${TITLES_ALIGNED}" "-DCONTAINS=>150<")

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
set(NO_PLOT_56 "${PLOT_T5}|${PLOT_T6}|No mapped reads.|Mismatches not available")

# report:na_notes. Unaligned input (FASTQ, FASTA, uBAM without @SQ): plots 1-4 only, no MAPQ
# table, no text of plots 5-6, no MD or identity notes.
plot_set_test(na_notes:fastq "${LONGREAD}/long.fq"
  "-DTITLES=${TITLES_UNALIGNED}" "-DNOT_TITLES=${TITLES_ABSENT_UNALIGNED}"
  "-DNOT_CONTAINS=No mapped reads.|Mismatches not available|Identity|% mapped|MAPQ >= 30"
  "-DSVG_MIN=4")
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
# (plot 5 keeps its svg), the notes name the reads that were left out
plot_set_test(na_notes:no_md "${SAMSTAT_FIXTURE_DIR}/hardclip_softclip.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=Mismatches not available: no MD tag.|Mismatch profile excludes 2 of 2 mapped reads|Insertions (% of aligned bases)|Deletions (% of aligned bases)")
# NM missing: the identity note
plot_set_test(na_notes:nonm "${LONGREAD}/nonm.sam"
  "-DTITLES=${TITLES_ALIGNED}"
  "-DCONTAINS=Identity unavailable for ")

# report:svg_a11y: every figure has role="img", <title> and <desc>
plot_set_test(svg_a11y:aligned "${SAMSTAT_FIXTURE_DIR}/paired_small.sam" "-DSVG_A11Y=1" "-DSVG_MIN=14")
plot_set_test(svg_a11y:ont_like "${LONGREAD}/ont_like.bam" "-DSVG_A11Y=1" "-DSVG_MIN=10")
plot_set_test(svg_a11y:fastq "${LONGREAD}/long.fq" "-DSVG_A11Y=1" "-DSVG_MIN=4")

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
  "-DCONTAINS=<h1>ont_like.bam</h1>|BAM &middot|long layout (median length ~1,422 bp|long when &gt;= 1 kb) &middot|SAMStat "
  "-DTITLES=<dt>Reads</dt>|<dt>Bases</dt>|<dt>N50</dt>|<dt>Mean Q</dt>|<dt>% mapped</dt>|<dt>% MAPQ &gt;= 30</dt>|<dt>% bases clipped</dt>|<dt>Median identity</dt>")
plot_set_test(header_line:aligned_short "${SAMSTAT_FIXTURE_DIR}/paired_small.sam"
  "-DCONTAINS=short layout (median length |<dt>Reads (R1+R2)</dt>|<th scope=\"col\" class=\"num\">Read 1</th>")
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
