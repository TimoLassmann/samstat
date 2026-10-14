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
