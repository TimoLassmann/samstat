# Reports are written to a temporary file .<name>.XXXXXX and renamed (src/report/outfile.c).
# Required: -DSAMSTAT=<binary> -DINPUT=<small file> -DOUTDIR=<dir> -DMODE=<mode>
#   MODE=symlink   every report name is a symlink (to a file, or dangling): samstat replaces
#                  the links with regular files, the link targets are unchanged / not created
#   MODE=nodir     no -d: the reports (and their temporary files) go next to the input
#   MODE=fail      a directory has the name of the JSON report: exit 1, the HTML report written
#                  before is complete, no MultiQC file is written
# Every mode: no temporary file (.*samstat*) is left in the output directory.

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
get_filename_component(name ${INPUT} NAME)
set(reports ${name}.samstat.html ${name}.samstat.json
  ${name}.samstat_general_mqc.json ${name}.samstat_mapq_mqc.json)

function(no_leftovers dir)
  file(GLOB left LIST_DIRECTORIES true "${dir}/.*samstat*")
  if(left)
    message(FATAL_ERROR "temporary files left behind: ${left}")
  endif()
endfunction()

function(run_samstat expect_rc)
  execute_process(
    TIMEOUT 120 COMMAND ${SAMSTAT} ${ARGN}
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)
  if(NOT "${rc}" STREQUAL "${expect_rc}")
    message(FATAL_ERROR "samstat ${ARGN}\nexit code ${rc}, expected ${expect_rc}\nstdout:\n${out}\nstderr:\n${err}")
  endif()
  set(log "${out}\n${err}" PARENT_SCOPE)
endfunction()

if(MODE STREQUAL "symlink")
  set(victim "${OUTDIR}/victim.txt")
  set(ghost "${OUTDIR}/ghost.txt")
  file(WRITE ${victim} "not a report\n")
  # the HTML and MultiQC names point to an existing file, the JSON name to a missing one
  foreach(r ${reports})
    if(r MATCHES "\\.samstat\\.json$")
      file(CREATE_LINK ${ghost} "${OUTDIR}/${r}" SYMBOLIC)
    else()
      file(CREATE_LINK ${victim} "${OUTDIR}/${r}" SYMBOLIC)
    endif()
  endforeach()
  run_samstat(0 --json --multiqc -d ${OUTDIR} ${INPUT})
  file(READ ${victim} v)
  if(NOT v STREQUAL "not a report\n")
    message(FATAL_ERROR "the symlink target was overwritten:\n${v}")
  endif()
  if(EXISTS ${ghost})
    message(FATAL_ERROR "a dangling symlink was followed: ${ghost} was created")
  endif()
  foreach(r ${reports})
    if(IS_SYMLINK "${OUTDIR}/${r}" OR NOT EXISTS "${OUTDIR}/${r}")
      message(FATAL_ERROR "${r} is not a regular file")
    endif()
  endforeach()
  file(READ "${OUTDIR}/${name}.samstat.html" html)
  if(NOT html MATCHES "</html>\n$")
    message(FATAL_ERROR "${name}.samstat.html is not a complete report")
  endif()
  no_leftovers(${OUTDIR})
elseif(MODE STREQUAL "nodir")
  file(MAKE_DIRECTORY ${OUTDIR}/in)
  file(COPY_FILE ${INPUT} ${OUTDIR}/in/${name})
  run_samstat(0 --json --multiqc ${OUTDIR}/in/${name})
  foreach(r ${reports})
    if(NOT EXISTS "${OUTDIR}/in/${r}")
      message(FATAL_ERROR "missing ${OUTDIR}/in/${r}")
    endif()
  endforeach()
  no_leftovers(${OUTDIR}/in)
  no_leftovers(${OUTDIR})
elseif(MODE STREQUAL "fail")
  file(MAKE_DIRECTORY "${OUTDIR}/${name}.samstat.json")
  run_samstat(1 --json --multiqc -d ${OUTDIR} ${INPUT})
  if(NOT log MATCHES "Failed to open output file")
    message(FATAL_ERROR "expected 'Failed to open output file' in the output:\n${log}")
  endif()
  file(READ "${OUTDIR}/${name}.samstat.html" html)
  if(NOT html MATCHES "</html>\n$")
    message(FATAL_ERROR "${name}.samstat.html is not a complete report")
  endif()
  file(GLOB mqc "${OUTDIR}/*_mqc.json")
  if(mqc)
    message(FATAL_ERROR "MultiQC files written after the JSON report failed: ${mqc}")
  endif()
  no_leftovers(${OUTDIR})
else()
  message(FATAL_ERROR "unknown MODE '${MODE}'")
endif()
