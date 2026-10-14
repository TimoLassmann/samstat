# SAMStat test suite.
#
# Test labels:
#   smoke      - samstat runs on every input in dev/ and writes a non-empty report
#
# Helper scripts live in tests/cmake/ and are run with `cmake -P` so the suite
# needs nothing beyond CMake and the samstat binary.

set(SAMSTAT_TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/dev)
set(SAMSTAT_TEST_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/out)
file(MAKE_DIRECTORY ${SAMSTAT_TEST_OUT_DIR})

# Smoke tests: one per input file in dev/.
file(GLOB SAMSTAT_SMOKE_INPUTS
  ${SAMSTAT_TEST_DATA_DIR}/*.sam
  ${SAMSTAT_TEST_DATA_DIR}/*.fa)

foreach(input ${SAMSTAT_SMOKE_INPUTS})
  get_filename_component(name ${input} NAME)
  set(outdir ${SAMSTAT_TEST_OUT_DIR}/smoke/${name})
  add_test(NAME smoke:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat>
      -DINPUT=${input}
      -DOUTDIR=${outdir}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_smoke.cmake)
  set_tests_properties(smoke:${name} PROPERTIES LABELS smoke)
endforeach()
