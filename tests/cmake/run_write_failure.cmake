# samstat must report a clear error (and write nothing) when the report cannot
# be written, instead of silently ignoring fopen/fwrite failures.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

# Running as root ignores directory permissions: nothing to test then.
execute_process(COMMAND id -u OUTPUT_VARIABLE uid OUTPUT_STRIP_TRAILING_WHITESPACE)
if(uid STREQUAL "0")
  message(STATUS "running as root, skipping")
  return()
endif()

file(CHMOD ${OUTDIR} PERMISSIONS OWNER_READ OWNER_EXECUTE)
execute_process(
  COMMAND ${SAMSTAT} -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
file(CHMOD ${OUTDIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)

# since 2.4 the output directory is checked before any record is read
string(FIND "${out}${err}" "use -d DIR" pos)
if(pos EQUAL -1)
  message(FATAL_ERROR "expected 'use -d DIR' in the output:\n${out}\n${err}")
endif()
file(GLOB written ${OUTDIR}/*)
if(written)
  message(FATAL_ERROR "files were written to a read-only directory: ${written}")
endif()
