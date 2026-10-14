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
