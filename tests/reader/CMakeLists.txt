# SAM/BAM reader tests (label: reader).
#
# reader:unit_parse_alignment  parse_alignment()/fix_orientation() on hand-built records
# reader:json:<fixture>        hand-computed values from tests/data (see tests/data/README.md)
# reader:error:<fixture>       unreadable / truncated input is rejected with a message
#
# Reads that overflow a heap buffer only fail visibly in a sanitizer build
# (cmake -DSAMSTAT_SANITIZE=ON), the value checks fail in every build.

set(READER_HELPERS ${PROJECT_SOURCE_DIR}/tests/cmake)
set(READER_OUT ${SAMSTAT_TEST_OUT_DIR}/reader)

# --- unit test -------------------------------------------------------------
add_executable(unit_parse_alignment
  unit_parse_alignment.c
  ${PROJECT_SOURCE_DIR}/src/htsinterface/htsglue.c
  ${PROJECT_SOURCE_DIR}/src/sambamparse/sam_bam_parse.c)
if(SAMSTAT_USE_SYSTEM_HTSLIB)
  target_link_libraries(unit_parse_alignment PRIVATE PkgConfig::HTSLIB)
else()
  add_dependencies(unit_parse_alignment htslib)
  target_include_directories(unit_parse_alignment PRIVATE ${HTSLIB_SRC_DIR})
  target_link_libraries(unit_parse_alignment PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)
endif()
if(Threads_FOUND)
  target_link_libraries(unit_parse_alignment PRIVATE Threads::Threads)
endif()
target_link_libraries(unit_parse_alignment PRIVATE tld-dev)
if(SAMSTAT_SANITIZE)
  target_compile_options(unit_parse_alignment PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(unit_parse_alignment PRIVATE -fsanitize=address,undefined)
endif()

add_test(NAME reader:unit_parse_alignment COMMAND unit_parse_alignment)
set_tests_properties(reader:unit_parse_alignment PROPERTIES LABELS reader)

# --- JSON value tests ------------------------------------------------------
# reader_json_test(<name> <fixture> <checks...>)
function(reader_json_test name input)
  # ';' inside an add_test argument has to be escaped to survive as one argument
  string(REPLACE ";" "\;" checks "${ARGN}")
  add_test(NAME reader:json:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${SAMSTAT_FIXTURE_DIR}/${input}
      -DOUTDIR=${READER_OUT}/${name} "-DCHECKS=${checks}"
      -P ${READER_HELPERS}/run_json.cmake)
  set_tests_properties(reader:json:${name} PROPERTIES LABELS reader)
endfunction()

# 1024 bases: exactly the size of the initial conversion buffer
reader_json_test(read_1024 read_1024.sam
  "summary.total_reads=1"
  "read1.length.mapq_ge30.1024=1")

# SEQ '*' with 10S40M: no length -10/40, no mismatches
reader_json_test(seq_star_softclip seq_star_softclip.sam
  "summary.total_reads=1"
  "summary.mapq_counts.read1.mapq_ge30=1"
  "read1.length.mapq_ge30#len=1"
  "read1.length.mapq_ge30.0=1")

# secondary (0x100) and supplementary (0x800) records are not reads
reader_json_test(supplementary supplementary.sam
  "summary.total_reads=1"
  "summary.secondary=1"
  "summary.supplementary=1"
  "summary.mapq_counts.read1.mapq_ge30=1"
  "summary.mapq_counts.read1.mapq_lt10=0"
  "read1.length.mapq_ge30.50=1")

# 3 primary reads, 2x flag 256 + 1x flag 2304 (secondary wins) + 1x flag 2048
reader_json_test(secondary secondary.sam
  "summary.total_reads=3"
  "summary.secondary=3"
  "summary.supplementary=1"
  "read1.n_reads=3")

# 5'S: the A->C mismatch is at aligned offset 9 on the forward read; the reverse
# read is reverse complemented (ref T, read G) so it is at offset 20-1-9 = 10.
reader_json_test(softclip5_md softclip5_md.sam
  "summary.total_reads=2"
  "read1.length.mapq_ge30.20=2"
  "read1.mismatches.mapq_ge30.C.9=1"
  "read1.mismatches.mapq_ge30.G.10=1"
  "read1.mismatches.mapq_ge30.A.9=0"
  "read1.mismatches.mapq_ge30.T.9=0")

# 3S: the read pointer must advance by the clip length, not by 1. The mismatch
# at aligned offset 2 is read base G (SEQ index 5), not A (SEQ index 3).
reader_json_test(softclip_shift softclip_shift.sam
  "summary.total_reads=1"
  "read1.length.mapq_ge30.20=1"
  "read1.mismatches.mapq_ge30.G.2=1"
  "read1.mismatches.mapq_ge30.A.2=0")

# 5H10S40M and 40M10S5H: only the soft clip is removed from the length (50 - 10)
reader_json_test(hardclip_softclip hardclip_softclip.sam
  "summary.total_reads=2"
  "read1.length.mapq_ge30.40=2")

# MD tags that disagree with the alignment: bounded by the alignment length
reader_json_test(md_inconsistent md_inconsistent.sam
  "summary.total_reads=3"
  "read1.length.mapq_ge30.20=3"
  "read1.mismatches.mapq_ge30.G.10=1"
  "read1.mismatches.mapq_ge30.G#len=11")
reader_json_test(md_8digit md_8digit.sam
  "summary.total_reads=1"
  "read1.length.mapq_ge30.20=1"
  "read1.mismatches.mapq_ge30.A#len=0")
reader_json_test(md_long_number md_long_number.sam
  "summary.total_reads=2"
  "read1.length.mapq_ge30.20=2"
  "read1.mismatches.mapq_ge30.G.10=1")

# 300M reads: the 256 byte slack of the alignment buffers is too small to hide an
# unbounded MD walk (MD digits 99999 / 300A99999 against 300 aligned bases)
reader_json_test(md_long_alignment md_long_alignment.sam
  "summary.total_reads=3"
  "read1.length.mapq_ge30.300=3"
  "read1.mismatches.mapq_ge30.A#len=0")

# zero byte file: no header can be read, still a valid empty input
reader_json_test(zero_bytes zero_bytes.sam
  "summary.total_reads=0")

# --- rejected inputs -------------------------------------------------------
# reader_error_test(<name> <fixture> <regex>)
function(reader_error_test name input regex)
  add_test(NAME reader:error:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${SAMSTAT_FIXTURE_DIR}/${input}
      -DOUTDIR=${READER_OUT}/err_${name} "-DEXPECT_REGEX=${regex}"
      -P ${READER_HELPERS}/run_expect_error.cmake)
  set_tests_properties(reader:error:${name} PROPERTIES LABELS reader)
endfunction()

reader_error_test(not_a_bam not_a_bam.bam "Cannot read the SAM/BAM header")
reader_error_test(truncated truncated.bam "truncated or corrupt")

# a 5000 bp read after three 20 bp reads: length histogram and quality rows grow with it
reader_json_test(long_after_short long_after_short.sam
  "summary.total_reads=4"
  "read1.length.mapq_ge30.20=3"
  "read1.length.mapq_ge30.5000=1"
  "read1.base_quality.mapq_ge30.40.19=4"
  "read1.base_quality.mapq_ge30.40.4999=1"
  "read1.base_quality.mapq_ge30.40#len=5000")

# --- reader:truncation_loop (WP2.1) --------------------------------------------------------
# A BAM cut off at arbitrary offsets must be rejected cleanly: exit code 1 or 2 (never 0,
# never a signal), a message, and no report left behind. The input is tests/data/longread/
# ont_like.bam (47 kB, two BGZF blocks plus the EOF marker). The cut offsets are chosen
# inside blocks; a cut exactly on a block boundary would be a valid shorter file.
# The cut copies are made with head(1) and deleted again. The script is generated into
# the build tree.
set(READER_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${READER_GEN})
file(WRITE ${READER_GEN}/run_truncation.cmake [=[
# -DSAMSTAT -DINPUT=<bam> -DOUTDIR -DCUTS=<offsets, '|' separated; negative = bytes before the end>
cmake_minimum_required(VERSION 3.19)
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
file(SIZE ${INPUT} total)
string(REPLACE "|" ";" cuts "${CUTS}")
set(failed "")
foreach(cut IN LISTS cuts)
  if(cut LESS 0)
    math(EXPR len "${total}${cut}")
  else()
    set(len ${cut})
  endif()
  set(dir ${OUTDIR}/cut${len})
  file(MAKE_DIRECTORY ${dir})
  set(bam ${dir}/cut.bam)
  execute_process(COMMAND head -c ${len} ${INPUT} OUTPUT_FILE ${bam} RESULT_VARIABLE hrc)
  if(NOT hrc EQUAL 0)
    message(FATAL_ERROR "head -c ${len} failed (${hrc})")
  endif()
  execute_process(COMMAND ${SAMSTAT} --json -d ${dir} ${bam}
    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc MATCHES "^[0-9]+$")
    string(APPEND failed "  cut at ${len}: crashed (${rc})\n")
  elseif(NOT (rc EQUAL 1 OR rc EQUAL 2))
    string(APPEND failed "  cut at ${len}: exit ${rc}, expected 1 or 2\n${out}${err}\n")
  elseif("${out}${err}" STREQUAL "")
    string(APPEND failed "  cut at ${len}: no message\n")
  endif()
  file(GLOB reports ${dir}/*.html ${dir}/*.json)
  if(reports)
    string(APPEND failed "  cut at ${len}: report left behind: ${reports}\n")
  endif()
  file(REMOVE ${bam})
endforeach()
if(failed)
  message(FATAL_ERROR "truncation check failed:\n${failed}")
endif()
]=])
add_test(NAME reader:truncation_loop
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DINPUT=${SAMSTAT_FIXTURE_DIR}/longread/ont_like.bam -DOUTDIR=${READER_OUT}/truncation_loop
    "-DCUTS=10|100|5000|20000|-100" -P ${READER_GEN}/run_truncation.cmake)
set_tests_properties(reader:truncation_loop PROPERTIES LABELS reader)
