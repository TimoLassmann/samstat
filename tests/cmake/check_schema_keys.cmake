# Run samstat --json on one input and compare the member names of the report with
# docs/samstat.schema.json: every 'required' member must be present and every
# member present must be listed in 'properties' (the schema has
# additionalProperties: false). Checked at the top level, in meta, summary,
# summary.mapq_counts and read1 (read2 when present).
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir> -DSCHEMA=<schema file>

cmake_minimum_required(VERSION 3.19)

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
execute_process(COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc}\n${out}\n${err}")
endif()
get_filename_component(name ${INPUT} NAME)
file(READ ${OUTDIR}/${name}.samstat.json doc)
file(READ ${SCHEMA} schema)

# member names of a JSON object, as a list
function(members obj out)
  string(JSON n LENGTH "${obj}")
  set(res "")
  if(n GREATER 0)
    math(EXPR last "${n}-1")
    foreach(i RANGE 0 ${last})
      string(JSON k MEMBER "${obj}" ${i})
      list(APPEND res ${k})
    endforeach()
  endif()
  set(${out} "${res}" PARENT_SCOPE)
endfunction()

# compare one object level against a schema node having 'properties'/'required'
function(check_level label schema_node doc_node)
  string(JSON props GET "${schema_node}" properties)
  members("${props}" prop_names)
  members("${doc_node}" doc_names)
  string(JSON nreq LENGTH "${schema_node}" required)
  if(nreq GREATER 0)
    math(EXPR lastreq "${nreq}-1")
    foreach(i RANGE 0 ${lastreq})
      string(JSON r GET "${schema_node}" required ${i})
      if(NOT r IN_LIST doc_names)
        message(SEND_ERROR "${INPUT}: ${label}: required member '${r}' missing from the report")
      endif()
    endforeach()
  endif()
  foreach(d IN LISTS doc_names)
    if(NOT d IN_LIST prop_names)
      message(SEND_ERROR "${INPUT}: ${label}: member '${d}' is not listed in the schema")
    endif()
  endforeach()
endfunction()

function(sub obj)
  string(JSON v GET "${obj}" ${ARGN})
  set(sub_result "${v}" PARENT_SCOPE)
endfunction()

check_level("top level" "${schema}" "${doc}")
string(JSON s_meta GET "${schema}" properties meta)
string(JSON d_meta GET "${doc}" meta)
check_level("meta" "${s_meta}" "${d_meta}")
string(JSON s_sum GET "${schema}" properties summary)
string(JSON d_sum GET "${doc}" summary)
check_level("summary" "${s_sum}" "${d_sum}")
string(JSON s_mc GET "${schema}" properties summary properties mapq_counts)
string(JSON d_mc GET "${doc}" summary mapq_counts)
check_level("summary.mapq_counts" "${s_mc}" "${d_mc}")
string(JSON s_read GET "${schema}" "$defs" read)
string(JSON d_r1 GET "${doc}" read1)
check_level("read1" "${s_read}" "${d_r1}")
string(JSON d_r2 ERROR_VARIABLE e2 GET "${doc}" read2)
if(e2 STREQUAL "NOTFOUND")
  check_level("read2" "${s_read}" "${d_r2}")
endif()
message(STATUS "schema keys of ${name} match")
