# Baseline golden tests (label: golden).
#
# For every dev/*.sam, dev/*.fa, tests/data/*.{sam,bam} and tests/data/longread/*.{sam,bam,fq}
# that ran without crashing when the goldens were captured, samstat --json is run and the report is compared
# with tests/golden/expected/<input>.json (meta ignored). The goldens capture
# CURRENT behaviour, including known bugs; when a fix intentionally changes a
# report, regenerate and review the diff:
#
#   cmake --build <build> --target update-goldens
#
# This is the only golden mechanism: json:golden:* tests (tests/json) compare against
# the same tests/golden/expected/<input>.json files.
#
# There is no size limit: the v2 collector keeps reports small for any read length.
#
# A golden test is registered only when tests/golden/expected/<input>.json exists. Goldens
# written before the "v2" object existed are compared without it (run_golden.cmake) until
# gate G3 regenerates them.

set(GOLDEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/expected)
set(GOLDEN_OUT ${SAMSTAT_TEST_OUT_DIR}/golden)

file(GLOB GOLDEN_CANDIDATES
  ${SAMSTAT_TEST_DATA_DIR}/*.sam ${SAMSTAT_TEST_DATA_DIR}/*.fa
  ${SAMSTAT_FIXTURE_DIR}/*.sam ${SAMSTAT_FIXTURE_DIR}/*.bam
  ${SAMSTAT_FIXTURE_DIR}/longread/*.sam ${SAMSTAT_FIXTURE_DIR}/longread/*.bam
  ${SAMSTAT_FIXTURE_DIR}/longread/*.fq)

set(update_commands)
foreach(input ${GOLDEN_CANDIDATES})
  get_filename_component(name ${input} NAME)
  set(cmd_args
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${input}"
    "-DOUTDIR=${GOLDEN_OUT}/${name}" "-DEXPECTED=${GOLDEN_DIR}/${name}.json")
  if(EXISTS "${GOLDEN_DIR}/${name}.json")
    add_test(NAME golden:${name}
      COMMAND ${CMAKE_COMMAND} ${cmd_args}
        -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_golden.cmake)
    set_tests_properties(golden:${name} PROPERTIES LABELS golden)
  endif()
  list(APPEND update_commands
    COMMAND ${CMAKE_COMMAND} ${cmd_args} -DUPDATE=ON
      -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_golden.cmake)
endforeach()

add_custom_target(update-goldens ${update_commands}
  COMMENT "Regenerating tests/golden/expected from the current samstat binary"
  VERBATIM)
add_dependencies(update-goldens samstat)
