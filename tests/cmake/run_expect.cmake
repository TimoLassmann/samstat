# Run samstat and check its exit code, output and the files it writes.
# Required: -DSAMSTAT=<binary> -DOUTDIR=<scratch dir, wiped first>
# Optional:
#   -DARGS=<list>                 command line; @OUTDIR@ and @GEN@ are substituted
#   -DEXPECT_RC=<n>               expected exit code (default 0)
#   -DEXPECT_STDERR_REGEX=<re>    must match the log output (stdout and stderr combined;
#                                 samstat writes all diagnostics to stderr)
#   -DEXPECT_NO_REGEX=<re>        must NOT match the log output
#   -DEXPECT_FILES=<list>         files (relative to OUTDIR) that must exist and be non-empty
#   -DEXPECT_NO_FILES=<list>      files (relative to OUTDIR) that must not exist
#   -DGEN_SAM=<n>                 first write OUTDIR/gen_<n>.sam with n reads, available as @GEN@
#   -DGEN_FASTQ=<n>               likewise OUTDIR/gen_<n>.fq with n reads, available as @GENFQ@
#   -DGEN_BLOB=<n>                OUTDIR/blob_<n>.bin: n bytes 'x', no newline, available as @GENBLOB@
#   -DJSON=<file rel. OUTDIR> -DCHECKS=<list>   forwarded to check_json.cmake
#   -DCOPY_FROM=<file> -DCOPY_AS=<name>   copy <file> to OUTDIR/<name> first, available as @COPY@
#   -DFILE_TEXT=<text> -DFILE_TEXT_IN=<file rel. OUTDIR>   the file must contain <text> literally
#   In ARGS, COPY_AS, EXPECT_FILES and FILE_TEXT: @Q@ is '"' and @BS@ is a backslash.
#   -DMKDIR=<list>                directories (rel. OUTDIR) created before the run, e.g. to make an output file name unwritable;
#                                 at most 100 characters and 4 levels each
#   -DNO_NAN=ON                   no 'nan' / 'inf' token in the files listed in EXPECT_FILES

cmake_minimum_required(VERSION 3.19)
get_filename_component(here "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

include(${CMAKE_CURRENT_LIST_DIR}/outdir_guard.cmake)
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

if(GEN_FASTQ)
  set(content "")
  set(chunk "")
  math(EXPR last "${GEN_FASTQ}-1")
  foreach(i RANGE 0 ${last})
    string(APPEND chunk "@q${i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")
    math(EXPR m "${i}%1000")
    if(m EQUAL 999)
      string(APPEND content "${chunk}")
      set(chunk "")
    endif()
  endforeach()
  string(APPEND content "${chunk}")
  set(GENFQ ${OUTDIR}/gen_${GEN_FASTQ}.fq)
  file(WRITE ${GENFQ} "${content}")
endif()

if(GEN_BLOB)
  string(REPEAT "x" ${GEN_BLOB} blob)
  set(GENBLOB ${OUTDIR}/blob_${GEN_BLOB}.bin)
  file(WRITE ${GENBLOB} "${blob}")
endif()

# characters that are awkward to pass through add_test
function(expand_tokens var)
  string(REPLACE "@Q@" "\"" v "${${var}}")
  string(REPLACE "@BS@" "\\" v "${v}")
  set(${var} "${v}" PARENT_SCOPE)
endfunction()

# Deep or very long paths below the build tree break backup tools (Time Machine) and
# `rm -rf` on some systems: MKDIR is for a few short names only (tests/README.md).
foreach(d IN LISTS MKDIR)
  string(LENGTH "${d}" dlen)
  string(REGEX MATCHALL "/" dslashes "${d}")
  list(LENGTH dslashes ddepth)
  if(dlen GREATER 100 OR ddepth GREATER 3)
    message(FATAL_ERROR "MKDIR '${d}' is too long or too deep for a test directory")
  endif()
  file(MAKE_DIRECTORY "${OUTDIR}/${d}")
endforeach()

if(COPY_FROM)
  expand_tokens(COPY_AS)
  file(COPY_FILE ${COPY_FROM} "${OUTDIR}/${COPY_AS}")
  set(COPY "${OUTDIR}/${COPY_AS}")
endif()

string(REPLACE "\;" ";" ARGS "${ARGS}")
set(run_args "")
foreach(a IN LISTS ARGS)
  expand_tokens(a)
  string(REPLACE "@COPY@" "${COPY}" a "${a}")
  string(REPLACE "@OUTDIR@" "${OUTDIR}" a "${a}")
  string(REPLACE "@GEN@" "${GEN}" a "${a}")
  string(REPLACE "@GENFQ@" "${GENFQ}" a "${a}")
  string(REPLACE "@GENBLOB@" "${GENBLOB}" a "${a}")
  list(APPEND run_args "${a}")
endforeach()

execute_process(
  TIMEOUT 120 COMMAND ${SAMSTAT} ${run_args}
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
  expand_tokens(f)
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
if(FILE_TEXT)
  expand_tokens(FILE_TEXT)
  expand_tokens(FILE_TEXT_IN)
  file(READ "${OUTDIR}/${FILE_TEXT_IN}" body)
  string(FIND "${body}" "${FILE_TEXT}" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "${FILE_TEXT_IN} does not contain '${FILE_TEXT}':\n${body}")
  endif()
endif()
string(REPLACE "\;" ";" EXPECT_NO_FILES "${EXPECT_NO_FILES}")
foreach(f IN LISTS EXPECT_NO_FILES)
  if(EXISTS ${OUTDIR}/${f})
    message(FATAL_ERROR "file ${f} must not exist")
  endif()
endforeach()

if(JSON)
  expand_tokens(JSON)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -DJSON=${OUTDIR}/${JSON} "-DCHECKS=${CHECKS}" -P ${here}/check_json.cmake
    RESULT_VARIABLE rc2 ERROR_VARIABLE err2 OUTPUT_VARIABLE out2)
  if(NOT rc2 EQUAL 0)
    message(FATAL_ERROR "${out2}${err2}")
  endif()
endif()
