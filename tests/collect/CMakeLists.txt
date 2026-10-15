# Statistics collection tests (label: collect).
#
# collect:qual:*        missing quality (QUAL '*') contributes nothing; qualities up to Q93 are kept
# collect:bounds:*      soft/hard clips, SEQ '*', end and normalised plots stay inside their arrays
# collect:counters      summary.secondary / summary.supplementary exist
# collect:plots:<file>  every dev/*.sam with --plotend --norm: exit 0 and valid JSON
# collect:mapq_sum:*    MAPQ class counts and length histograms add up to the primary reads
# collect:axis_rule:*   layout rule and both position axes
# collect:unit          unit_collect: the v2 collector (merge, tables, position bins, per-record
#                       decisions, 128-bit mean Q, summary values of the long-read fixtures)
#
# Expected values come from the hand-computed fixtures in tests/data/README.md.

# --- unit test: the collector compiled into its own executable ---------------
# collect.c includes collect/convert_tables.h (committed to the source tree, next to collect.c).
# gen/ holds the check scripts written further down.
set(COLLECT_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${COLLECT_GEN_DIR})
add_executable(unit_collect
  unit_collect.c
  ${PROJECT_SOURCE_DIR}/src/collect/collect.c
  ${PROJECT_SOURCE_DIR}/src/htsinterface/htsglue.c
  ${PROJECT_SOURCE_DIR}/src/sambamparse/sam_bam_parse.c)
target_compile_definitions(unit_collect PRIVATE "COLLECT_TEST_DATA_DIR=\"${SAMSTAT_FIXTURE_DIR}\"")
target_link_libraries(unit_collect PRIVATE htslib::hts)
if(Threads_FOUND)
  target_link_libraries(unit_collect PRIVATE Threads::Threads)
endif()
target_link_libraries(unit_collect PRIVATE module_plot tldlite)
if(UNIX)
  target_link_libraries(unit_collect PRIVATE m)
endif()
if(SAMSTAT_SANITIZE)
  target_compile_options(unit_collect PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(unit_collect PRIVATE -fsanitize=address,undefined)
endif()
add_test(NAME collect:unit COMMAND unit_collect)
set_tests_properties(collect:unit PROPERTIES LABELS collect)

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
  read1.length.reads.mapq_ge30.30=1 read1.length.reads.mapq_10_29.30=1 read1.length.reads.mapq_lt10.30=1
  "read1.absolute.qual_bases.mapq_ge30#len=0"
  "read1.absolute.qual_bases.mapq_10_29#len=0"
  "read1.absolute.qual_bases.mapq_lt10#len=0"
  "read1.relative.qual_error_sum.mapq_ge30#len=0"
  "read1.read_mean_q.mapq_ge30#len=0"
  "read1.read_mean_q.mapq_10_29#len=0"
  "read1.read_mean_q.mapq_lt10#len=0"
  summary.mean_q=)
collect_test(qual:missing_default ${SAMSTAT_FIXTURE_DIR}/no_qual.sam
  CHECKS ${NO_QUAL_COMMON})
collect_test(qual:missing_plotend_norm ${SAMSTAT_FIXTURE_DIR}/no_qual.sam
  ARGS --plotend --norm
  CHECKS ${NO_QUAL_COMMON}
  HTML_ABSENT "41\\.0000")

# high_qual.sam: hq1 = Q40 x5, Q41 x5, Q42 x5, Q60 x5, Q93 x5; hq2 = Q60 x25.
# 2.4 keeps every value up to Q93 (the 2.3 cap at Q41 is gone). E[q] = round(10^(-q/10) * 2^30):
# E[40]=107374 E[41]=85290 E[42]=67749 E[60]=1074 E[93]=1. Error sums per position:
# 0..4: E40+E60 = 108448; 5..9: E41+E60 = 86364; 10..14: E42+E60 = 68823;
# 15..19: 2 * E60 = 2148; 20..24: E93+E60 = 1075. Per-read mean Q: hq2 is Q60;
# hq1 has mean error 1307440 / 25 = 52297, the largest q with E[q] >= 52297 is 43.
collect_test(qual:clamp_high ${SAMSTAT_FIXTURE_DIR}/high_qual.sam
  ARGS --plotend --norm
  CHECKS
    summary.total_reads=2
    "read1.absolute.qual_bases.mapq_ge30#len=25"
    read1.absolute.qual_bases.mapq_ge30.0=2
    read1.absolute.qual_bases.mapq_ge30.24=2
    read1.absolute.qual_error_sum.mapq_ge30.0=108448
    read1.absolute.qual_error_sum.mapq_ge30.4=108448
    read1.absolute.qual_error_sum.mapq_ge30.5=86364
    read1.absolute.qual_error_sum.mapq_ge30.10=68823
    read1.absolute.qual_error_sum.mapq_ge30.15=2148
    read1.absolute.qual_error_sum.mapq_ge30.20=1075
    read1.absolute.qual_error_sum.mapq_ge30.24=1075
    "read1.read_mean_q.mapq_ge30#len=61"
    read1.read_mean_q.mapq_ge30.43=1
    read1.read_mean_q.mapq_ge30.60=1)

# --- index bounds ----------------------------------------------------------
# seq_star_softclip.sam: SEQ '*' with CIGAR 10S40M. The stored length is 0, not
# 0 - 10 (out-of-bounds write before the fix).
collect_test(bounds:seq_star_softclip ${SAMSTAT_FIXTURE_DIR}/seq_star_softclip.sam
  ARGS --plotend --norm
  CHECKS
    summary.total_reads=1
    read1.length.reads.mapq_ge30.0=1
    "read1.absolute.qual_bases.mapq_ge30#len=0")

# hardclip_softclip.sam: 50 bp reads with 10 soft + 5 hard clip, 40 bp aligned: the position
# axes cover the aligned part (40 positions, Q40 on both reads).
collect_test(bounds:hardclip_softclip ${SAMSTAT_FIXTURE_DIR}/hardclip_softclip.sam
  ARGS --plotend --norm
  CHECKS
    summary.total_reads=2
    "read1.absolute.qual_bases.mapq_ge30#len=40"
    read1.absolute.qual_bases.mapq_ge30.0=2
    read1.absolute.qual_bases.mapq_ge30.39=2)

collect_test(bounds:softclip5_md ${SAMSTAT_FIXTURE_DIR}/softclip5_md.sam
  ARGS --plotend --norm
  CHECKS summary.total_reads=2)

# 4 reads from 2 pairs: read2 plots and R2 end / norm plots are exercised.
collect_test(bounds:paired_plotend_norm ${SAMSTAT_FIXTURE_DIR}/paired_small.sam
  ARGS --plotend --norm
  CHECKS
    summary.total_reads=4
    summary.read1_reads=2
    summary.read2_reads=2
    read2.n_reads=2)

collect_test(bounds:small_bam_plotend_norm ${SAMSTAT_FIXTURE_DIR}/small.bam
  ARGS --plotend --norm
  CHECKS summary.total_reads=20)

# --- counters --------------------------------------------------------------
# The reader fills the counters; the collection must at least report them (as 0
# here) in summary.
collect_test(counters ${SAMSTAT_FIXTURE_DIR}/small.sam
  CHECKS summary.secondary=0 summary.supplementary=0 summary.total_reads=20)

# --- every dev/*.sam with --plotend --norm ---------------------------------
file(GLOB COLLECT_DEV_SAM ${SAMSTAT_TEST_DATA_DIR}/*.sam)
foreach(input ${COLLECT_DEV_SAM})
  get_filename_component(name ${input} NAME)
  collect_test(plots:${name} ${input} ARGS --plotend --norm)
endforeach()

# --- format 2 summary and axes -----------------------------------------------
# collect:mapq_sum:<input>  summary.mapq_counts add up to the primary reads, and so do the
#                           per-class length histograms of the read blocks
# collect:axis_rule:<len>   reads of <len> bp: layout short (median < 1000) or long, both
#                           position axes present
file(MAKE_DIRECTORY ${COLLECT_GEN_DIR}/scripts)
file(WRITE ${COLLECT_GEN_DIR}/scripts/check_v2_mapq_sum.cmake [=[
# -DSAMSTAT -DINPUT -DOUTDIR
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\n${out}\n${err}")
endif()
get_filename_component(name ${INPUT} NAME)
file(READ ${OUTDIR}/${name}.samstat.json doc)
# string(JSON) re-parses its input on every call: take the parts apart first
string(JSON sum_obj GET "${doc}" summary)
string(JSON total GET "${sum_obj}" total_reads)
string(JSON r2 GET "${sum_obj}" read2_reads)
set(mates read1)
if(r2 GREATER 0)
  list(APPEND mates read2)
endif()
set(sum 0)
set(len_sum 0)
foreach(m IN LISTS mates)
  foreach(c mapq_ge30 mapq_10_29 mapq_lt10 unmapped)
    string(JSON v GET "${sum_obj}" mapq_counts ${m} ${c})
    math(EXPR sum "${sum}+${v}")
    string(JSON arr GET "${doc}" ${m} length reads ${c})
    string(JSON n LENGTH "${arr}")
    if(n GREATER 0)
      math(EXPR last "${n}-1")
      foreach(i RANGE ${last})
        string(JSON v GET "${arr}" ${i})
        math(EXPR len_sum "${len_sum}+${v}")
      endforeach()
    endif()
  endforeach()
endforeach()
if(NOT sum EQUAL total)
  message(FATAL_ERROR "${name}: MAPQ class counts add up to ${sum}, total_reads is ${total}")
endif()
if(NOT len_sum EQUAL total)
  message(FATAL_ERROR "${name}: length histograms add up to ${len_sum}, total_reads is ${total}")
endif()
]=])
foreach(input ${SAMSTAT_FIXTURE_DIR}/small.sam ${SAMSTAT_FIXTURE_DIR}/paired_small.sam
    ${SAMSTAT_FIXTURE_DIR}/longread/ont_like.sam ${SAMSTAT_FIXTURE_DIR}/longread/single_len.sam
    ${SAMSTAT_FIXTURE_DIR}/longread/supp_chain.sam ${SAMSTAT_FIXTURE_DIR}/longread/sq_nomap.sam)
  get_filename_component(name ${input} NAME)
  add_test(NAME collect:mapq_sum:${name}
    COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${input}
      -DOUTDIR=${COLLECT_OUT}/mapq_sum_${name}
      -P ${COLLECT_GEN_DIR}/scripts/check_v2_mapq_sum.cmake)
  set_tests_properties(collect:mapq_sum:${name} PROPERTIES LABELS collect)
endforeach()

file(WRITE ${COLLECT_GEN_DIR}/scripts/check_axis_rule.cmake [=[
# -DSAMSTAT -DOUTDIR -DLEN=<read length> -DLAYOUT=short|long
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
math(EXPR quarter "${LEN}/4+1")
string(REPEAT "ACGT" ${quarter} seq)
string(SUBSTRING "${seq}" 0 ${LEN} seq)
string(REPEAT "I" ${LEN} qual)
set(sam ${OUTDIR}/axis.sam)
file(WRITE ${sam} "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n")
foreach(i 1 2 3)
  file(APPEND ${sam} "r${i}\t0\tchr1\t${i}\t60\t${LEN}M\t*\t0\t0\t${seq}\t${qual}\tNM:i:0\tMD:Z:${LEN}\n")
endforeach()
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${sam}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc}\n${out}\n${err}")
endif()
file(READ ${OUTDIR}/axis.sam.samstat.json doc)
string(JSON layout GET "${doc}" summary layout)
if(NOT layout STREQUAL LAYOUT)
  message(FATAL_ERROR "reads of ${LEN} bp: layout '${layout}', expected '${LAYOUT}'")
endif()
# both axes are filled whatever the layout; the first read base of every read is an A
# (the relative axis has 1% bins, so bin 0 holds the first bases of several positions)
foreach(axis absolute relative)
  string(JSON n LENGTH "${doc}" read1 ${axis} composition mapq_ge30 A)
  if(n LESS 1)
    message(FATAL_ERROR "reads of ${LEN} bp: read1.${axis} is empty")
  endif()
  string(JSON a0 GET "${doc}" read1 ${axis} composition mapq_ge30 A 0)
  if(axis STREQUAL "absolute" AND NOT a0 EQUAL 3)
    message(FATAL_ERROR "reads of ${LEN} bp: read1.${axis} bin 0 has ${a0} A, expected 3")
  elseif(a0 LESS 3)
    message(FATAL_ERROR "reads of ${LEN} bp: read1.${axis} bin 0 has ${a0} A, expected at least 3")
  endif()
endforeach()
]=])
foreach(spec "150;short" "999;short" "1000;long" "10000;long")
  list(GET spec 0 len)
  list(GET spec 1 layout)
  add_test(NAME collect:axis_rule:${len}
    COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${COLLECT_OUT}/axis_rule_${len}
      -DLEN=${len} -DLAYOUT=${layout} -P ${COLLECT_GEN_DIR}/scripts/check_axis_rule.cmake)
  set_tests_properties(collect:axis_rule:${len} PROPERTIES LABELS collect)
endforeach()

# --- batch memory bounds -----------------------------------------------------------
# collect:batch_bytes:threads  records with large aux tags, long reads (20000, 30000 and
#                              32767 bp) and a CIGAR of 9000 operations, cut into batches by
#                              the record byte limit (2 * SAMSTAT_BATCH_BASES = 40000 bytes:
#                              four short reads with a 12 KB tag, before 8 records or 20000
#                              bases) or the base limit (the long reads), with slot buffers
#                              trimmed on capacity: -t 1, 2, 3, 4 and 8 give identical results
# collect:batch_bytes:counts   the same input without the test-only limits: every read counted
# The input (about 300 KB) is written here, at configure time, into the build tree.
set(BATCH_BYTES_SAM ${COLLECT_GEN_DIR}/batch_bytes.sam)
string(REPEAT "A" 12000 big_tag)
string(REPEAT "ACGTTGCA" 19 short_seq)
string(SUBSTRING "${short_seq}" 0 150 short_seq)
string(REPEAT "I" 150 short_qual)
set(batch_bytes_text "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n")
set(batch_bytes_long 20000 30000 32767 20000)
set(i 0)
foreach(k RANGE 1 8)
  string(APPEND batch_bytes_text "s${k}\t0\tchr1\t${k}\t60\t150M\t*\t0\t0\t${short_seq}\t${short_qual}\tNM:i:0\tMD:Z:150\tXB:Z:${big_tag}\n")
  if(k EQUAL 4)
    # 4500 x 1M1I: 9000 bases, 9000 CIGAR operations
    string(REPEAT "1M1I" 4500 cig)
    string(REPEAT "ACGTTGCA" 1125 seq)
    string(REPEAT "I" 9000 qual)
    string(APPEND batch_bytes_text "c1\t0\tchr1\t100\t60\t${cig}\t*\t0\t0\t${seq}\t${qual}\tNM:i:4500\tMD:Z:4500\n")
  endif()
  if(k LESS 5)
    list(GET batch_bytes_long ${i} len)
    math(EXPR i "${i}+1")
    math(EXPR reps "${len}/8+1")
    string(REPEAT "ACGTTGCA" ${reps} seq)
    string(SUBSTRING "${seq}" 0 ${len} seq)
    string(REPEAT "I" ${len} qual)
    string(APPEND batch_bytes_text "l${i}\t0\tchr1\t${k}\t60\t${len}M\t*\t0\t0\t${seq}\t${qual}\tNM:i:0\tMD:Z:${len}\n")
  endif()
endforeach()
file(WRITE ${BATCH_BYTES_SAM} "${batch_bytes_text}")
unset(batch_bytes_text)
unset(big_tag)
unset(seq)
unset(qual)
unset(cig)
add_test(NAME collect:batch_bytes:threads
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${BATCH_BYTES_SAM}
    -DOUTDIR=${COLLECT_OUT}/batch_bytes_threads -DMODE=one_file
    -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_threads.cmake)
set_tests_properties(collect:batch_bytes:threads PROPERTIES LABELS collect)
collect_test(batch_bytes:counts ${BATCH_BYTES_SAM}
  CHECKS summary.total_reads=13)
