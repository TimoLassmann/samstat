# One read of 1 Mb must give a small JSON: all arrays are bounded, independent of read length.
# Writes a single-read SAM (1,000,000 bp, CIGAR 1000000M, NM and MD present) into OUTDIR
# with string(REPEAT), runs samstat --json and checks (standard-plots.md section 7):
#   JSON file < 1 MB; every absolute-axis array <= 1000 elements; every relative-axis array
#   <= 100; identity arrays <= 1001; read_mean_q <= 94.
# Required: -DSAMSTAT=<binary> -DOUTDIR=<dir>

cmake_minimum_required(VERSION 3.19)

include(${CMAKE_CURRENT_LIST_DIR}/outdir_guard.cmake)
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

set(len 1000000)
math(EXPR quarter "${len}/4")
string(REPEAT "ACGT" ${quarter} seq)
string(REPEAT "I" ${len} qual)
set(sam ${OUTDIR}/ultralong.sam)
file(WRITE ${sam}
  "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000000\n"
  "ultra\t0\tchr1\t1\t60\t${len}M\t*\t0\t0\t${seq}\t${qual}\tNM:i:0\tMD:Z:${len}\n")
file(SIZE ${sam} sam_size)
if(sam_size LESS ${len})
  message(FATAL_ERROR "generated SAM is only ${sam_size} bytes")
endif()

execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${sam}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on the 1 Mb read\nstdout:\n${out}\nstderr:\n${err}")
endif()
set(report ${OUTDIR}/ultralong.sam.samstat.json)
if(NOT EXISTS ${report})
  message(FATAL_ERROR "No JSON report written: ${report}")
endif()

file(SIZE ${report} json_size)
if(NOT json_size LESS 1000000)
  message(FATAL_ERROR "JSON report is ${json_size} bytes, expected < 1000000")
endif()
file(READ ${report} doc)

# walk(<limit> <json path...>): every array below the path has at most <limit> elements
set(failed "")
function(walk limit)
  set(path ${ARGN})
  string(JSON type ERROR_VARIABLE err TYPE "${doc}" ${path})
  if(NOT err STREQUAL "NOTFOUND")
    string(APPEND failed "  missing: ${path}\n")
  elseif(type STREQUAL "OBJECT")
    string(JSON n LENGTH "${doc}" ${path})
    if(n GREATER 0)
      math(EXPR last "${n}-1")
      foreach(i RANGE ${last})
        string(JSON member MEMBER "${doc}" ${path} ${i})
        walk(${limit} ${path} ${member})
      endforeach()
    endif()
  elseif(type STREQUAL "ARRAY")
    string(JSON n LENGTH "${doc}" ${path})
    if(n GREATER ${limit})
      list(JOIN path "." p)
      string(APPEND failed "  ${p}: ${n} elements, limit ${limit}\n")
    endif()
  endif()
  set(failed "${failed}" PARENT_SCOPE)
endfunction()

foreach(spec "read1;absolute;1000" "read1;relative;100" "read1;identity;1001" "read1;read_mean_q;94")
  list(GET spec 0 a)
  list(GET spec 1 b)
  list(GET spec 2 limit)
  set(start ${a} ${b})
  # the block must exist and must not be empty (a test that walks nothing proves nothing)
  string(JSON n ERROR_VARIABLE err LENGTH "${doc}" ${start})
  if(NOT err STREQUAL "NOTFOUND" OR n EQUAL 0)
    list(JOIN start "." p)
    string(APPEND failed "  ${p} is missing or empty (${err})\n")
    continue()
  endif()
  walk(${limit} ${start})
endforeach()

if(failed)
  message(FATAL_ERROR "ultralong JSON check failed (${report}, ${json_size} bytes):\n${failed}")
endif()
message(STATUS "ultralong: JSON ${json_size} bytes, all array bounds hold")
