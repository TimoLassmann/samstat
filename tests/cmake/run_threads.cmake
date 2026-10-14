# The thread count must not change the results: 8 copies of one input, processed with
# -t 1 and -t 4, give the same JSON (meta ignored). Also runs once without -d and with
# --verbose: the report lands next to the input.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR}/in)
set(copies "")
foreach(i RANGE 1 8)
  file(COPY_FILE ${INPUT} ${OUTDIR}/in/s${i}.sam)
  list(APPEND copies ${OUTDIR}/in/s${i}.sam)
endforeach()

foreach(t 1 4)
  file(MAKE_DIRECTORY ${OUTDIR}/t${t})
  execute_process(COMMAND ${SAMSTAT} -t ${t} --json -d ${OUTDIR}/t${t} ${copies}
    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "samstat -t ${t} exited with ${rc}\n${out}\n${err}")
  endif()
  foreach(i RANGE 1 8)
    foreach(ext samstat.html samstat.json)
      if(NOT EXISTS ${OUTDIR}/t${t}/s${i}.sam.${ext})
        message(FATAL_ERROR "-t ${t}: s${i}.sam.${ext} missing")
      endif()
    endforeach()
  endforeach()
endforeach()

foreach(i RANGE 1 8)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -DACTUAL=${OUTDIR}/t4/s${i}.sam.samstat.json
      -DEXPECTED=${OUTDIR}/t1/s${i}.sam.samstat.json -P ${here}/compare_golden.cmake
    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "s${i}.sam: -t 1 and -t 4 differ\n${out}\n${err}")
  endif()
endforeach()

# no -d: report next to the input; --verbose is accepted
execute_process(COMMAND ${SAMSTAT} --verbose ${OUTDIR}/in/s1.sam
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0 OR NOT EXISTS ${OUTDIR}/in/s1.sam.samstat.html)
  message(FATAL_ERROR "run without -d failed (rc ${rc})\n${out}\n${err}")
endif()
