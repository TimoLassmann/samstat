# Run samstat on one input and validate the JSON report.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>
# Optional: -DNO_JSON=ON    run WITHOUT --json and assert no .json file is written
#           -DCHECKS=<list> forwarded to check_json.cmake
#           -DGOLDEN=<file> compared with compare_golden.cmake (meta ignored);
#                           add -DUPDATE=ON to regenerate it

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

set(args -d ${OUTDIR})
if(NOT NO_JSON)
  list(APPEND args --json)
endif()

execute_process(
  COMMAND ${SAMSTAT} ${args} ${INPUT}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\nstdout:\n${out}\nstderr:\n${err}")
endif()

if(NOT EXISTS ${OUTDIR}/${name}.samstat.html)
  message(FATAL_ERROR "No HTML report written for ${INPUT}")
endif()

if(NO_JSON)
  file(GLOB stray ${OUTDIR}/*.json)
  if(stray)
    message(FATAL_ERROR "JSON written without --json: ${stray}")
  endif()
  return()
endif()

if(NOT EXISTS ${report})
  message(FATAL_ERROR "No JSON report written: ${report}")
endif()

# string(JSON) re-parses the whole document on every call. Reports of long-read
# inputs (e.g. dev/ce#large_seq.sam: ~100 MB) would take minutes, so for those only
# the framing is checked.
file(SIZE ${report} report_size)
if(report_size GREATER 4000000)
  message(STATUS "${report} is ${report_size} bytes: checking framing only")
  file(READ ${report} head LIMIT 512)
  if(NOT head MATCHES "^{\n\"format_version\": 1,\n\"meta\": {\"samstat_version\": ")
    message(FATAL_ERROR "${report}: unexpected start: ${head}")
  endif()
  math(EXPR tail_off "${report_size}-4")
  file(READ ${report} tail OFFSET ${tail_off})
  if(NOT tail STREQUAL "}\n}\n")
    message(FATAL_ERROR "${report}: unexpected end: '${tail}'")
  endif()
  return()
endif()

file(READ ${report} doc)
string(JSON version ERROR_VARIABLE err GET "${doc}" format_version)
if(NOT err STREQUAL "NOTFOUND")
  message(FATAL_ERROR "${report} is not valid JSON or lacks format_version: ${err}")
endif()
if(NOT version STREQUAL "1")
  message(FATAL_ERROR "Unexpected format_version '${version}'")
endif()
foreach(key meta summary mapq_bins read1)
  string(JSON ignored ERROR_VARIABLE err GET "${doc}" ${key})
  if(NOT err STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${report}: missing top-level key '${key}': ${err}")
  endif()
endforeach()
foreach(key total_reads read1_reads read2_reads paired properly_paired is_aligned is_partial_report mapq_counts)
  string(JSON ignored ERROR_VARIABLE err GET "${doc}" summary ${key})
  if(NOT err STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${report}: missing summary.${key}: ${err}")
  endif()
endforeach()
foreach(key n_reads length base_composition base_quality mismatches insertions deletions)
  string(JSON ignored ERROR_VARIABLE err GET "${doc}" read1 ${key})
  if(NOT err STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${report}: missing read1.${key}: ${err}")
  endif()
endforeach()
string(JSON n_bins LENGTH "${doc}" mapq_bins)
if(NOT n_bins EQUAL 4)
  message(FATAL_ERROR "${report}: expected 4 mapq_bins, got ${n_bins}")
endif()
# read2 is present exactly when there are read 2 records
string(JSON n_r2 GET "${doc}" summary read2_reads)
string(JSON ignored ERROR_VARIABLE err GET "${doc}" read2)
if(n_r2 GREATER 0 AND NOT err STREQUAL "NOTFOUND")
  message(FATAL_ERROR "${report}: read2_reads=${n_r2} but no read2 section")
elseif(n_r2 EQUAL 0 AND err STREQUAL "NOTFOUND")
  message(FATAL_ERROR "${report}: read2 section present but read2_reads=0")
endif()
# per-bin read counts must add up to total_reads
string(JSON total GET "${doc}" summary total_reads)
set(sum 0)
foreach(r read1 read2)
  string(JSON ignored ERROR_VARIABLE err GET "${doc}" summary mapq_counts ${r})
  if(err STREQUAL "NOTFOUND")
    foreach(b mapq_ge30 mapq_10_29 mapq_lt10 unmapped)
      string(JSON v GET "${doc}" summary mapq_counts ${r} ${b})
      math(EXPR sum "${sum}+${v}")
    endforeach()
  endif()
endforeach()
if(NOT sum EQUAL total)
  message(FATAL_ERROR "${report}: mapq_counts sum ${sum} != total_reads ${total}")
endif()

string(REPLACE "\\;" ";" CHECKS "${CHECKS}")
if(CHECKS)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -DJSON=${report} "-DCHECKS=${CHECKS}" -P ${here}/check_json.cmake
    RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_VARIABLE out)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${out}${err}")
  endif()
endif()

if(GOLDEN)
  if(NOT UPDATE)
    golden_without_v2(${report} ${GOLDEN} compare_report)
  else()
    set(compare_report ${report})
  endif()
  execute_process(
    COMMAND ${CMAKE_COMMAND} -DACTUAL=${compare_report} -DEXPECTED=${GOLDEN} -DUPDATE=${UPDATE} -P ${here}/compare_golden.cmake
    RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_VARIABLE out)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${out}${err}")
  endif()
endif()
