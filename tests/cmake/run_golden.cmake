# Run samstat --json on one input and compare the report with a golden file.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir> -DEXPECTED=<file>
# Optional: -DUPDATE=ON  overwrite EXPECTED with the new report instead of comparing.
#                        Inputs on which samstat fails or whose report exceeds
#                        MAX_BYTES (default 4000000) are skipped silently.

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)
if(NOT MAX_BYTES)
  set(MAX_BYTES 4000000)
endif()

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
get_filename_component(name ${INPUT} NAME)
set(report ${OUTDIR}/${name}.samstat.json)

execute_process(
  COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)

if(UPDATE)
  if(NOT rc EQUAL 0 OR NOT EXISTS ${report})
    message(STATUS "skip ${name}: samstat exit ${rc} / no JSON report")
    return()
  endif()
  file(SIZE ${report} sz)
  if(sz GREATER MAX_BYTES)
    message(STATUS "skip ${name}: report is ${sz} bytes")
    return()
  endif()
else()
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\n${out}\n${err}")
  endif()
  if(NOT EXISTS ${report})
    message(FATAL_ERROR "No JSON report written: ${report}")
  endif()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -DACTUAL=${report} -DEXPECTED=${EXPECTED} -DUPDATE=${UPDATE}
    -P ${here}/compare_golden.cmake
  RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_VARIABLE out)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${out}${err}")
endif()
