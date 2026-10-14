# Run samstat (optionally with --multiqc) on one input and check which MultiQC
# custom content files appear.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>
# Optional: -DMULTIQC=ON (pass --multiqc and require all 4 files)
#           -DMULTIQC=OFF (do not pass it and require that no *_mqc.json exists)
# The output directory is NOT cleaned when -DKEEP=ON, so several inputs can share it.

if(NOT KEEP)
  file(REMOVE_RECURSE ${OUTDIR})
endif()
file(MAKE_DIRECTORY ${OUTDIR})

set(args -d ${OUTDIR})
if(MULTIQC)
  list(APPEND args --multiqc)
endif()

execute_process(
  COMMAND ${SAMSTAT} ${args} ${INPUT}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\nstdout:\n${out}\nstderr:\n${err}")
endif()

get_filename_component(name ${INPUT} NAME)
file(GLOB found ${OUTDIR}/${name}.samstat_*_mqc.json)
if(MULTIQC)
  foreach(section general mapq length basequal)
    if(NOT EXISTS ${OUTDIR}/${name}.samstat_${section}_mqc.json)
      message(FATAL_ERROR "Missing ${OUTDIR}/${name}.samstat_${section}_mqc.json")
    endif()
  endforeach()
else()
  if(found)
    message(FATAL_ERROR "_mqc files written without --multiqc: ${found}")
  endif()
  if(NOT EXISTS ${OUTDIR}/${name}.samstat.html)
    message(FATAL_ERROR "No HTML report written")
  endif()
endif()
