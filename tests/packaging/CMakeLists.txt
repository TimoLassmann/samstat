add_test(NAME packaging:metadata
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPROJECT_VERSION=${PROJECT_VERSION}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/check_packaging.cmake)
set_tests_properties(packaging:metadata PROPERTIES LABELS packaging)

# Fixture budget: tests/data must stay at or below 1 MB (see docs/plan-v2.4, T1).
set(FIXTURE_BUDGET_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/fixture_budget.cmake)
file(WRITE ${FIXTURE_BUDGET_SCRIPT} [=[
file(GLOB_RECURSE files LIST_DIRECTORIES false "${DIR}/*")
set(total 0)
foreach(f IN LISTS files)
  file(SIZE "${f}" sz)
  math(EXPR total "${total}+${sz}")
endforeach()
if(total GREATER ${MAX_BYTES})
  message(FATAL_ERROR "fixture_budget: ${DIR} holds ${total} bytes, budget is ${MAX_BYTES}")
endif()
message(STATUS "fixture_budget: ${total} of ${MAX_BYTES} bytes used in ${DIR}")
]=])
add_test(NAME packaging:fixture_budget
  COMMAND ${CMAKE_COMMAND} -DDIR=${SAMSTAT_FIXTURE_DIR} -DMAX_BYTES=1000000
    -P ${FIXTURE_BUDGET_SCRIPT})
set_tests_properties(packaging:fixture_budget PROPERTIES LABELS packaging)
