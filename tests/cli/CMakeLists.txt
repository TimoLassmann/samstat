# Command line, threading and exit code tests (label: cli).
#
# Exit codes: 0 success (also for empty inputs, with a warning), 1 at least one
# input failed (the other inputs are still processed), 2 usage error.
# All tests use tests/cmake/run_expect.cmake.

set(CLI_HELPERS ${PROJECT_SOURCE_DIR}/tests/cmake)
set(CLI_OUT ${SAMSTAT_TEST_OUT_DIR}/cli)

# cli_test(<name> ARGS <args...> [RC n] [REGEX re] [NO_REGEX re] [FILES f...]
#          [NO_FILES f...] [GEN_SAM n] [JSON file CHECKS c...] [NO_NAN])
# Arguments may use @OUTDIR@ (scratch dir of the test) and @GEN@ (generated SAM).
function(cli_test name)
  cmake_parse_arguments(T "NO_NAN" "RC;REGEX;NO_REGEX;GEN_SAM;JSON" "ARGS;FILES;NO_FILES;CHECKS" ${ARGN})
  if(NOT DEFINED T_RC)
    set(T_RC 0)
  endif()
  string(REPLACE ";" "\;" args "${T_ARGS}")
  string(REPLACE ";" "\;" files "${T_FILES}")
  string(REPLACE ";" "\;" nofiles "${T_NO_FILES}")
  string(REPLACE ";" "\;" checks "${T_CHECKS}")
  set(cmd ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${CLI_OUT}/${name}
    "-DARGS=${args}" -DEXPECT_RC=${T_RC})
  if(T_REGEX)
    list(APPEND cmd "-DEXPECT_STDERR_REGEX=${T_REGEX}")
  endif()
  if(T_NO_REGEX)
    list(APPEND cmd "-DEXPECT_NO_REGEX=${T_NO_REGEX}")
  endif()
  if(T_FILES)
    list(APPEND cmd "-DEXPECT_FILES=${files}")
  endif()
  if(T_NO_FILES)
    list(APPEND cmd "-DEXPECT_NO_FILES=${nofiles}")
  endif()
  if(T_GEN_SAM)
    list(APPEND cmd -DGEN_SAM=${T_GEN_SAM})
  endif()
  if(T_JSON)
    list(APPEND cmd "-DJSON=${T_JSON}" "-DCHECKS=${checks}")
  endif()
  if(T_NO_NAN)
    list(APPEND cmd -DNO_NAN=ON)
  endif()
  list(APPEND cmd -P ${CLI_HELPERS}/run_expect.cmake)
  add_test(NAME cli:${name} COMMAND ${cmd})
  set_tests_properties(cli:${name} PROPERTIES LABELS cli)
endfunction()

set(FX ${SAMSTAT_FIXTURE_DIR})
set(D @OUTDIR@)

# --- worker error handling / exit codes of failing inputs ---------------------
# 1 = input failed, other inputs still processed
foreach(threads 1 4)
  cli_test(mixed_good_bad_t${threads}
    ARGS -t ${threads} -d ${D} ${FX}/small.sam ${FX}/does_not_exist.sam
    RC 1 REGEX "does_not_exist\\.sam"
    FILES small.sam.samstat.html)
  cli_test(mixed_bad_good_t${threads}
    ARGS -t ${threads} -d ${D} ${FX}/does_not_exist.sam ${FX}/small.sam
    RC 1 REGEX "does_not_exist\\.sam"
    FILES small.sam.samstat.html)
endforeach()
cli_test(nonexistent_file
  ARGS -d ${D} ${FX}/does_not_exist.bam
  RC 1 REGEX "does_not_exist\\.bam")
cli_test(all_good_rc0
  ARGS -t 2 -d ${D} ${FX}/small.sam ${FX}/small.bam
  FILES small.sam.samstat.html small.bam.samstat.html)
