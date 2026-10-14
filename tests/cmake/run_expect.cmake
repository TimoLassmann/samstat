# Run samstat and check its exit code, output and the files it writes.
# Required: -DSAMSTAT=<binary> -DOUTDIR=<scratch dir, wiped first>
# Optional:
#   -DARGS=<list>                 command line; @OUTDIR@ and @GEN@ are substituted
#   -DEXPECT_RC=<n>               expected exit code (default 0)
#   -DEXPECT_STDERR_REGEX=<re>    must match the log output (stdout and stderr combined,
#                                 samstat writes warnings to stdout and errors to stderr)
#   -DEXPECT_NO_REGEX=<re>        must NOT match the log output
#   -DEXPECT_FILES=<list>         files (relative to OUTDIR) that must exist and be non-empty
#   -DEXPECT_NO_FILES=<list>      files (relative to OUTDIR) that must not exist
#   -DGEN_SAM=<n>                 first write OUTDIR/gen_<n>.sam with n reads, available as @GEN@
#   -DJSON=<file rel. OUTDIR> -DCHECKS=<list>   forwarded to check_json.cmake
#   -DNO_NAN=ON                   no 'nan' / 'inf' token in the files listed in EXPECT_FILES

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

if(NOT DEFINED EXPECT_RC)
  set(EXPECT_RC 0)
endif()

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

if(GEN_SAM)
  set(gen ${OUTDIR}/gen_${GEN_SAM}.sam)
  set(content "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n")
  set(chunk "")
  math(EXPR last "${GEN_SAM}-1")
  foreach(i RANGE 0 ${last})
    math(EXPR pos "${i}%90000+1")
    string(APPEND chunk "r${i}\t0\tchr1\t${pos}\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n")
    math(EXPR m "${i}%1000")
    if(m EQUAL 999)
      string(APPEND content "${chunk}")
      set(chunk "")
    endif()
  endforeach()
  string(APPEND content "${chunk}")
  file(WRITE ${gen} "${content}")
  set(GEN ${gen})
endif()

string(REPLACE "\;" ";" ARGS "${ARGS}")
set(run_args "")
foreach(a IN LISTS ARGS)
  string(REPLACE "@OUTDIR@" "${OUTDIR}" a "${a}")
  string(REPLACE "@GEN@" "${GEN}" a "${a}")
  list(APPEND run_args "${a}")
endforeach()

execute_process(
  COMMAND ${SAMSTAT} ${run_args}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
set(log "${out}\n${err}")

if(NOT "${rc}" STREQUAL "${EXPECT_RC}")
  message(FATAL_ERROR "samstat ${run_args}\nexit code ${rc}, expected ${EXPECT_RC}\nstdout:\n${out}\nstderr:\n${err}")
endif()
if(EXPECT_STDERR_REGEX AND NOT log MATCHES "${EXPECT_STDERR_REGEX}")
  message(FATAL_ERROR "log output does not match '${EXPECT_STDERR_REGEX}'\nstdout:\n${out}\nstderr:\n${err}")
endif()
if(EXPECT_NO_REGEX AND log MATCHES "${EXPECT_NO_REGEX}")
  message(FATAL_ERROR "log output unexpectedly matches '${EXPECT_NO_REGEX}'\nstdout:\n${out}\nstderr:\n${err}")
endif()

string(REPLACE "\;" ";" EXPECT_FILES "${EXPECT_FILES}")
foreach(f IN LISTS EXPECT_FILES)
  if(NOT EXISTS ${OUTDIR}/${f})
    message(FATAL_ERROR "expected file ${f} was not written\nstdout:\n${out}\nstderr:\n${err}")
  endif()
  file(SIZE ${OUTDIR}/${f} sz)
  if(sz EQUAL 0)
    message(FATAL_ERROR "${f} is empty")
  endif()
  if(NO_NAN AND sz LESS 4000000)
    file(READ ${OUTDIR}/${f} body)
    if(body MATCHES "[^A-Za-z0-9_]-?(nan|inf)[^A-Za-z0-9_]")
      message(FATAL_ERROR "${f} contains a nan/inf token: ${CMAKE_MATCH_0}")
    endif()
  endif()
endforeach()
string(REPLACE "\;" ";" EXPECT_NO_FILES "${EXPECT_NO_FILES}")
foreach(f IN LISTS EXPECT_NO_FILES)
  if(EXISTS ${OUTDIR}/${f})
    message(FATAL_ERROR "file ${f} must not exist")
  endif()
endforeach()

if(JSON)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -DJSON=${OUTDIR}/${JSON} "-DCHECKS=${CHECKS}" -P ${here}/check_json.cmake
    RESULT_VARIABLE rc2 ERROR_VARIABLE err2 OUTPUT_VARIABLE out2)
  if(NOT rc2 EQUAL 0)
    message(FATAL_ERROR "${out2}${err2}")
  endif()
endif()
