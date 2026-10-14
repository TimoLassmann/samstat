# Run samstat on one input and check that it exits 0 and writes a non-empty report.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

execute_process(
  COMMAND ${SAMSTAT} -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)

if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\nstdout:\n${out}\nstderr:\n${err}")
endif()

get_filename_component(name ${INPUT} NAME)
set(report ${OUTDIR}/${name}.samstat.html)
if(NOT EXISTS ${report})
  message(FATAL_ERROR "No report written: ${report}")
endif()

file(SIZE ${report} size)
if(size EQUAL 0)
  message(FATAL_ERROR "Empty report: ${report}")
endif()
