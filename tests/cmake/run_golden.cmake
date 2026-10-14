# Run samstat --json on one input and compare the report with a golden file.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir> -DEXPECTED=<file>
# Optional: -DUPDATE=ON  overwrite EXPECTED with the new report instead of comparing.
#                        Inputs on which samstat fails are skipped silently.

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

# Transitional (WP3.1b to gate G3): the report carries a "v2" object that goldens captured
# before G3 do not have. Compare without it when the golden has none; once the golden has "v2"
# (after regeneration) it is compared in full.
function(golden_without_v2 report expected out)
  file(READ ${expected} exp)
  string(JSON ignored ERROR_VARIABLE err GET "${exp}" v2)
  set(${out} ${report} PARENT_SCOPE)
  if(err STREQUAL "NOTFOUND")
    return()
  endif()
  file(READ ${report} act)
  string(JSON act ERROR_VARIABLE err REMOVE "${act}" v2)
  if(NOT err STREQUAL "NOTFOUND")
    return()
  endif()
  file(WRITE ${report}.nov2 "${act}")
  set(${out} ${report}.nov2 PARENT_SCOPE)
endfunction()

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
else()
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\n${out}\n${err}")
  endif()
  if(NOT EXISTS ${report})
    message(FATAL_ERROR "No JSON report written: ${report}")
  endif()
endif()

if(NOT UPDATE)
  golden_without_v2(${report} ${EXPECTED} compare_report)
else()
  set(compare_report ${report})
endif()
execute_process(
  COMMAND ${CMAKE_COMMAND} -DACTUAL=${compare_report} -DEXPECTED=${EXPECTED} -DUPDATE=${UPDATE}
    -P ${here}/compare_golden.cmake
  RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_VARIABLE out)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${out}${err}")
endif()
