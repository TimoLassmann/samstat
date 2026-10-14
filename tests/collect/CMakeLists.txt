# Statistics collection tests (label: collect).
#
# collect:qual:*        missing quality (QUAL '*') contributes nothing; qualities are binned 0..41
# collect:bounds:*      soft/hard clips, SEQ '*', end and normalised plots stay inside their arrays
# collect:counters      summary.secondary / summary.supplementary exist
# collect:plots:<file>  every dev/*.sam with --plotend --norm: exit 0 and valid JSON
#
# Expected values come from the hand-computed fixtures in tests/data/README.md.

set(COLLECT_HELPER ${PROJECT_SOURCE_DIR}/tests/cmake/run_collect.cmake)
set(COLLECT_OUT ${SAMSTAT_TEST_OUT_DIR}/collect)

# collect_test(<name> <input> [ARGS <opts>...] [CHECKS <checks>...]
#              [HTML_ABSENT <regex>] [HTML_PRESENT <regex>])
function(collect_test name input)
  cmake_parse_arguments(T "" "HTML_ABSENT;HTML_PRESENT" "ARGS;CHECKS" ${ARGN})
  # ';' inside an add_test argument has to be escaped to survive as one argument
  string(REPLACE ";" "\\;" args "${T_ARGS}")
  string(REPLACE ";" "\\;" checks "${T_CHECKS}")
  add_test(NAME collect:${name}
    COMMAND ${CMAKE_COMMAND}
      "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${input}"
      "-DOUTDIR=${COLLECT_OUT}/${name}"
      "-DARGS=${args}" "-DCHECKS=${checks}"
      "-DHTML_ABSENT=${T_HTML_ABSENT}" "-DHTML_PRESENT=${T_HTML_PRESENT}"
      -P ${COLLECT_HELPER})
  set_tests_properties(collect:${name} PROPERTIES LABELS collect)
endfunction()

# --- missing quality -------------------------------------------------------
# no_qual.sam: 3 reads of 30 bp, QUAL '*'. Reads without quality must not appear
# in any quality plot. Before the fix --plotend / --norm recorded them as Q41.
set(NO_QUAL_COMMON
  summary.total_reads=3
  read1.length.mapq_ge30.30=1 read1.length.mapq_10_29.30=1 read1.length.mapq_lt10.30=1
  "read1.base_quality.mapq_ge30#len=0"
  "read1.base_quality.mapq_10_29#len=0"
  "read1.base_quality.mapq_lt10#len=0")
collect_test(qual:missing_default ${SAMSTAT_FIXTURE_DIR}/no_qual.sam
  CHECKS ${NO_QUAL_COMMON})
collect_test(qual:missing_plotend_norm ${SAMSTAT_FIXTURE_DIR}/no_qual.sam
  ARGS --plotend --norm
  CHECKS ${NO_QUAL_COMMON}
  HTML_ABSENT "41\\.0000")

# high_qual.sam: hq1 = Q40 x5, Q41 x5, Q42 x5, Q60 x5, Q93 x5; hq2 = Q60 x25.
# Everything above Q41 is counted in bin 41; nothing lands outside 0..41.
# Position 0..4: Q40 once (hq1), Q41 once (hq2). Position 5..24: Q41 twice.
collect_test(qual:clamp_high ${SAMSTAT_FIXTURE_DIR}/high_qual.sam
  ARGS --plotend --norm
  CHECKS
    summary.total_reads=2
    "read1.base_quality.mapq_ge30#len=42"
    read1.base_quality.mapq_ge30.40.0=1
    read1.base_quality.mapq_ge30.40.4=1
    read1.base_quality.mapq_ge30.40.5=0
    read1.base_quality.mapq_ge30.41.0=1
    read1.base_quality.mapq_ge30.41.4=1
    read1.base_quality.mapq_ge30.41.5=2
    read1.base_quality.mapq_ge30.41.24=2)

