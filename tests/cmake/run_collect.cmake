# Run samstat with extra options on one input, then check the JSON and HTML reports.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>
# Optional: -DARGS=<list>          extra samstat options, e.g. "--plotend;--norm"
#           -DCHECKS=<list>        forwarded to check_json.cmake
#           -DHTML_ABSENT=<regex>  the HTML report must NOT match
#           -DHTML_PRESENT=<regex> the HTML report must match
# --json is always added. The JSON must be valid; reports over 4 MB (long reads)
# only get a framing check because string(JSON) is slow on them.

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

file(REMOVE_RECURSE "${OUTDIR}")
file(MAKE_DIRECTORY "${OUTDIR}")

string(REPLACE "\\;" ";" ARGS "${ARGS}")
get_filename_component(name "${INPUT}" NAME)
set(json "${OUTDIR}/${name}.samstat.json")
set(html "${OUTDIR}/${name}.samstat.html")

execute_process(
  TIMEOUT 120 COMMAND ${SAMSTAT} ${ARGS} --json -d "${OUTDIR}" "${INPUT}"
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat ${ARGS} exited with ${rc} on ${INPUT}\nstdout:\n${out}\nstderr:\n${err}")
endif()
foreach(f "${json}" "${html}")
  if(NOT EXISTS "${f}")
    message(FATAL_ERROR "Report not written: ${f}")
  endif()
endforeach()

file(SIZE "${json}" json_size)
if(json_size GREATER 4000000)
  file(READ "${json}" head LIMIT 64)
  math(EXPR tail_off "${json_size}-4")
  file(READ "${json}" tail OFFSET ${tail_off})
  if(NOT head MATCHES "^{\n\"format_version\": 2," OR NOT tail STREQUAL "}\n}\n")
    message(FATAL_ERROR "${json}: bad framing")
  endif()
else()
  file(READ "${json}" doc)
  string(JSON version ERROR_VARIABLE jerr GET "${doc}" format_version)
  if(NOT jerr STREQUAL "NOTFOUND" OR NOT version STREQUAL "2")
    message(FATAL_ERROR "${json}: invalid JSON or format_version: ${jerr}")
  endif()
  if(CHECKS)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -DJSON=${json} "-DCHECKS=${CHECKS}" -P ${here}/check_json.cmake
      RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_VARIABLE out)
    if(NOT rc EQUAL 0)
      message(FATAL_ERROR "${out}${err}")
    endif()
  endif()
endif()

if(HTML_ABSENT OR HTML_PRESENT)
  file(READ "${html}" page)
  if(HTML_ABSENT AND page MATCHES "${HTML_ABSENT}")
    message(FATAL_ERROR "${html} must not match '${HTML_ABSENT}' (found '${CMAKE_MATCH_0}')")
  endif()
  if(HTML_PRESENT AND NOT page MATCHES "${HTML_PRESENT}")
    message(FATAL_ERROR "${html} does not match '${HTML_PRESENT}'")
  endif()
endif()
