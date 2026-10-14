# The thread count must not change the results: 8 copies of one input, processed with
# -t 1 and -t 4, give the same JSON (meta ignored). Also runs once without -d and with
# --verbose: the report lands next to the input.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>
# Optional: -DMODE=one_file (see below); the default is the many-inputs test above.

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

# -DMODE=one_file: ONE input, many threads (intra-file threading, WP4.1).
# Runs -t 1, 3 and 8 on INPUT with the test-only batch limits SAMSTAT_BATCH_RECORDS=8 and
# SAMSTAT_BATCH_BASES=20000, so that even a small file is cut into several batches, and
# requires byte-identical JSON. --verbose must report 'workers: N' with N > 1 for -t 3
# and -t 8. Optional -DREPLICATE=<n>: the records of a SAM INPUT are repeated n times, to
# reach at least 5 batches (a 4-record fixture cannot fill them).
if(MODE STREQUAL "one_file")
  file(REMOVE_RECURSE ${OUTDIR})
  file(MAKE_DIRECTORY ${OUTDIR})
  get_filename_component(iname ${INPUT} NAME)
  set(run_input ${INPUT})
  if(REPLICATE)
    file(STRINGS ${INPUT} lines)
    set(text "")
    set(records "")
    foreach(l IN LISTS lines)
      if(l MATCHES "^@")
        string(APPEND text "${l}\n")
      else()
        string(APPEND records "${l}\n")
      endif()
    endforeach()
    set(run_input ${OUTDIR}/${iname})
    file(WRITE ${run_input} "${text}")
    foreach(i RANGE 1 ${REPLICATE})
      file(APPEND ${run_input} "${records}")
    endforeach()
  endif()
  foreach(t 1 3 8)
    file(MAKE_DIRECTORY ${OUTDIR}/t${t})
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E env SAMSTAT_BATCH_RECORDS=8 SAMSTAT_BATCH_BASES=20000
        ${SAMSTAT} -t ${t} --verbose --json -d ${OUTDIR}/t${t} ${run_input}
      RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
      message(FATAL_ERROR "samstat -t ${t} exited with ${rc} on ${iname}\n${out}\n${err}")
    endif()
    if(NOT EXISTS ${OUTDIR}/t${t}/${iname}.samstat.json)
      message(FATAL_ERROR "-t ${t}: ${iname}.samstat.json missing\n${out}\n${err}")
    endif()
    if(t GREATER 1)
      if(NOT "${out}\n${err}" MATCHES "workers: ([0-9]+)")
        message(FATAL_ERROR "-t ${t} --verbose printed no 'workers: N' line\n${out}\n${err}")
      endif()
      if(NOT CMAKE_MATCH_1 GREATER 1)
        message(FATAL_ERROR "-t ${t}: workers: ${CMAKE_MATCH_1}, expected more than 1 on one file\n${out}\n${err}")
      endif()
    endif()
  endforeach()
  file(READ ${OUTDIR}/t1/${iname}.samstat.json j1)
  foreach(t 3 8)
    file(READ ${OUTDIR}/t${t}/${iname}.samstat.json jt)
    if(NOT j1 STREQUAL jt)
      message(FATAL_ERROR "${iname}: JSON differs between -t 1 and -t ${t}")
    endif()
  endforeach()
  if(REPLICATE)
    file(REMOVE ${run_input})
  endif()
  return()
endif()

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
