# Run samstat on one input and check the generated HTML report.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>
# Optional (lists are separated with '|' because ';' does not survive add_test):
#   -DAS_NAME=<name>        copy INPUT to OUTDIR/<name> first and run on the copy
#   -DCONTAINS=a|b          every item must occur literally in the report
#   -DNOT_CONTAINS=a|b      no item may occur literally in the report
#   -DMATCHES=re|re         every regular expression must match
#   -DNOT_MATCHES=re|re     no regular expression may match

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

set(run_input ${INPUT})
if(AS_NAME)
  set(run_input ${OUTDIR}/${AS_NAME})
  file(COPY_FILE ${INPUT} ${run_input})
endif()
get_filename_component(name ${run_input} NAME)

execute_process(
  TIMEOUT 120 COMMAND ${SAMSTAT} -d ${OUTDIR} ${run_input}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on ${run_input}\nstdout:\n${out}\nstderr:\n${err}")
endif()

set(report ${OUTDIR}/${name}.samstat.html)
if(NOT EXISTS ${report})
  message(FATAL_ERROR "No report written: ${report}\nstderr:\n${err}")
endif()
file(READ ${report} html)

set(failed "")
if(CONTAINS)
  string(REPLACE "|" ";" items "${CONTAINS}")
  foreach(item IN LISTS items)
    string(FIND "${html}" "${item}" pos)
    if(pos EQUAL -1)
      string(APPEND failed "  missing: ${item}\n")
    endif()
  endforeach()
endif()
if(NOT_CONTAINS)
  string(REPLACE "|" ";" items "${NOT_CONTAINS}")
  foreach(item IN LISTS items)
    string(FIND "${html}" "${item}" pos)
    if(NOT pos EQUAL -1)
      string(APPEND failed "  unexpected: ${item}\n")
    endif()
  endforeach()
endif()
if(MATCHES)
  string(REPLACE "|" ";" items "${MATCHES}")
  foreach(item IN LISTS items)
    if(NOT "${html}" MATCHES "${item}")
      string(APPEND failed "  no match for regex: ${item}\n")
    endif()
  endforeach()
endif()
if(NOT_MATCHES)
  string(REPLACE "|" ";" items "${NOT_MATCHES}")
  foreach(item IN LISTS items)
    if("${html}" MATCHES "${item}")
      string(APPEND failed "  unexpected regex match: ${item}\n")
    endif()
  endforeach()
endif()

if(failed)
  message(FATAL_ERROR "HTML check failed for ${report}:\n${failed}")
endif()
