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

# --- usage errors exit 2, help / version exit 0 --------------------------------
cli_test(no_args ARGS "" RC 2 REGEX "Usage")
cli_test(help_short ARGS -h RC 0 REGEX "Usage")
cli_test(help_long ARGS -help RC 0 REGEX "Usage")
cli_test(version_short ARGS -v RC 0 REGEX "[0-9]+\\.[0-9]+")
cli_test(version_long ARGS -version RC 0 REGEX "[0-9]+\\.[0-9]+")
cli_test(only_options ARGS -t 2 RC 2 REGEX "Usage")
cli_test(unknown_option ARGS --bogus ${FX}/small.sam RC 2 NO_FILES small.sam.samstat.html)
cli_test(missing_option_value ARGS ${FX}/small.sam -t RC 2)
cli_test(duplicate_input
  ARGS -d ${D} ${FX}/small.sam ${FX}/small.sam
  RC 2 REGEX "Duplicated input" NO_FILES small.sam.samstat.html)
cli_test(outdir_missing
  ARGS -d ${D}/nope ${FX}/small.sam RC 2 REGEX "does not exist")

# --- file type detection --------------------------------------------------------
cli_test(ext_unknown_only ARGS -d ${D} ${FX}/unknown.xyz RC 2 REGEX "unknown\\.xyz")
# an unrecognised extension is an error even when other inputs are valid
cli_test(ext_unknown_with_good
  ARGS -d ${D} ${FX}/small.sam ${FX}/unknown.xyz
  RC 2 REGEX "unknown\\.xyz" NO_FILES small.sam.samstat.html)
cli_test(ext_fq_gz ARGS -d ${D} ${FX}/reads.fq.gz FILES reads.fq.gz.samstat.html)
cli_test(ext_fa_gz ARGS -d ${D} ${FX}/reads.fa.gz FILES reads.fa.gz.samstat.html)
cli_test(ext_upper_case ARGS -d ${D} ${FX}/READS.FQ FILES READS.FQ.samstat.html)
cli_test(ext_not_a_suffix_only
  ARGS -d ${D} ${FX}/small.sam.txt RC 2 REGEX "small\\.sam\\.txt")

# upper/mixed case and the .gz variants of every fasta/fastq extension, via copies
set(EXT_COPY_DIR ${CLI_OUT}/ext_copies)
file(MAKE_DIRECTORY ${EXT_COPY_DIR})
foreach(ext SAM Sam BAM Bam)
  if(ext MATCHES "^[Ss][Aa][Mm]$")
    set(src ${FX}/small.sam)
  else()
    set(src ${FX}/small.bam)
  endif()
  configure_file(${src} ${EXT_COPY_DIR}/case.${ext} COPYONLY)
  cli_test(ext_case_${ext} ARGS -d ${D} ${EXT_COPY_DIR}/case.${ext} FILES case.${ext}.samstat.html)
endforeach()
foreach(ext fa fasta fna fq fastq FA FASTQ)
  configure_file(${FX}/READS.FQ ${EXT_COPY_DIR}/plain.${ext} COPYONLY)
  cli_test(ext_plain_${ext} ARGS -d ${D} ${EXT_COPY_DIR}/plain.${ext} FILES plain.${ext}.samstat.html)
endforeach()
foreach(ext fa.gz fasta.gz fna.gz fq.gz fastq.gz FQ.GZ Fa.Gz)
  configure_file(${FX}/reads.fq.gz ${EXT_COPY_DIR}/gz.${ext} COPYONLY)
  cli_test(ext_gz_${ext} ARGS -d ${D} ${EXT_COPY_DIR}/gz.${ext} FILES gz.${ext}.samstat.html)
endforeach()
