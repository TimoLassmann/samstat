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

report_html_test(plotly_pinned "${REPORT_SAM}"
  "-DNOT_CONTAINS=plotly-latest"
  "-DCONTAINS=https://cdn.plot.ly/plotly-1.58.5.min.js"
  "-DMATCHES=plotly-[0-9]+\\.[0-9]+\\.[0-9]+\\.min\\.js")
report_html_test(no_unused_d3 "${REPORT_SAM}"
  "-DNOT_CONTAINS=d3js.org")

# 'lines' is not a Plotly trace type; line traces must be scatter + mode lines.
foreach(input "xx#minimal.sam" "xx#pair.sam" "ce#1.sam" "ce#2.sam")
  report_html_test(trace_type:${input} "${SAMSTAT_TEST_DATA_DIR}/${input}"
    "-DNOT_CONTAINS=type:'lines'|type: 'lines'|type:\"lines\"|type: \"lines\"")
endforeach()
report_html_test(trace_type_lines_mode "${SAMSTAT_TEST_DATA_DIR}/ce#1.sam"
  "-DCONTAINS=type: 'scatter',\nmode: 'lines'")

# HTML hygiene
report_html_test(html_lang_viewport "${REPORT_SAM}"
  "-DCONTAINS=<html lang=\"en\">|<meta charset=\"utf-8\">|<meta name=\"viewport\""
  "-DNOT_CONTAINS=100%%")
# A file name with '<' and '&' must be escaped in <title> and <h1> (copy made at test time).
report_html_test(escape_filename "${REPORT_SAM}"
  "-DAS_NAME=a<b&c.sam"
  "-DCONTAINS=<h1>a&lt;b&amp;c.sam</h1>|a&lt;b&amp;c.sam\n</title>"
  "-DNOT_CONTAINS=a<b&c.sam</h1>|<title>\na<b")

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
# WP2.1: pending report tests of the 2.4 standard report. All DISABLED until WP4.2a
# (report page) enables them. Specification: docs/plan-v2.4/standard-plots.md sections 1, 2
# and 5.
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
  set_tests_properties(report:${name} PROPERTIES LABELS report DISABLED TRUE)   # enabled by WP4.2a
endfunction()

# report:plot_set (T5): the titles in order, nothing else; no NaN|Infinity|undefined
plot_set_test(plot_set:short_bam "${SAMSTAT_FIXTURE_DIR}/small.bam" "-DTITLES=${TITLES_ALIGNED}")
plot_set_test(plot_set:ont_like_bam "${SAMSTAT_FIXTURE_DIR}/longread/ont_like.bam" "-DTITLES=${TITLES_ALIGNED}")
plot_set_test(plot_set:long_fq "${SAMSTAT_FIXTURE_DIR}/longread/long.fq"
  "-DTITLES=${TITLES_UNALIGNED}" "-DNOT_TITLES=${TITLES_ABSENT_UNALIGNED}")
# the short layout x axis ends at the longest read (150 bp for single_len.sam), not at 1000.
# The exact tick markup is WP4.2a's: adjust the pattern there if the axis labels differ.
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
  PROPERTIES LABELS report DISABLED TRUE)   # enabled by WP4.2a
