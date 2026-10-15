# The table header collect.c includes is committed (src/collect/convert_tables.h); this runs the
# generator it came from and fails if the two differ, so that they cannot drift apart.
# Required: -DGEN=<converttablegen> -DOUT=<file in the build tree> -DREF=<committed header>
file(REMOVE "${OUT}")
execute_process(COMMAND "${GEN}" "${OUT}" RESULT_VARIABLE rc ERROR_VARIABLE err TIMEOUT 120)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "converttablegen failed (${rc}):\n${err}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${OUT}" "${REF}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${REF} differs from the output of converttablegen (${OUT}); "
    "if convert_tables.c changed on purpose, copy ${OUT} over ${REF}")
endif()
