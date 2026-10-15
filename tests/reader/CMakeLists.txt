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
target_link_libraries(unit_parse_alignment PRIVATE htslib::hts)
if(Threads_FOUND)
  target_link_libraries(unit_parse_alignment PRIVATE Threads::Threads)
endif()
target_link_libraries(unit_parse_alignment PRIVATE tldlite)
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

# 1024 bases: exactly the size of the initial conversion buffer. Storage length bins are exact
# below 1000 bp and then follow g(1000) = 1000, 1010, 1021, ...: 1024 is in the bin of edge 1021,
# which is bin 1002.
reader_json_test(read_1024 read_1024.sam
  "summary.total_reads=1"
  "summary.min_length=1024"
  "summary.max_length=1024"
  "read1.length.reads.mapq_ge30#len=1003"
  "read1.length.reads.mapq_ge30.1002=1"
  "read1.length.bases.mapq_ge30.1002=1024"
  "read1.absolute.qual_bases.mapq_ge30#len=1000")

# SEQ '*' with 10S40M: no length -10/40, no mismatches
reader_json_test(seq_star_softclip seq_star_softclip.sam
  "summary.total_reads=1"
  "summary.mapq_counts.read1.mapq_ge30=1"
  "read1.length.reads.mapq_ge30#len=1"
  "read1.length.reads.mapq_ge30.0=1"
  "read1.absolute.qual_bases.mapq_ge30#len=0"
  "read1.absolute.mismatches.mapq_ge30#len=0")

# secondary (0x100) and supplementary (0x800) records are not reads
reader_json_test(supplementary supplementary.sam
  "summary.total_reads=1"
  "summary.secondary=1"
  "summary.supplementary=1"
  "summary.mapq_counts.read1.mapq_ge30=1"
  "summary.mapq_counts.read1.mapq_lt10=0"
  "read1.length.reads.mapq_ge30.50=1")

# 3 primary reads, 2x flag 256 + 1x flag 2304 (secondary wins) + 1x flag 2048
reader_json_test(secondary secondary.sam
  "summary.total_reads=3"
  "summary.secondary=3"
  "summary.supplementary=1"
  "read1.n_reads=3")

# 5'S (5S20M, 25 bp): the A->C mismatch is at aligned offset 9 on the forward read; the reverse
# read is reverse complemented (ref T, read G) so it is at offset 20-1-9 = 10. The length is the
# full read length (25), the position axes cover the aligned part (20 positions, both reads have MD).
reader_json_test(softclip5_md softclip5_md.sam
  "summary.total_reads=2"
  "read1.length.reads.mapq_ge30.25=2"
  "read1.absolute.aligned_md.mapq_ge30#len=20"
  "read1.absolute.aligned_md.mapq_ge30.9=2"
  "read1.absolute.mismatches.mapq_ge30#len=11"
  "read1.absolute.mismatches.mapq_ge30.9=1"
  "read1.absolute.mismatches.mapq_ge30.10=1"
  "read1.absolute.mismatches.mapq_ge30.8=0"
  "read1.relative.mismatches.mapq_ge30.45=1"
  "read1.relative.mismatches.mapq_ge30.50=1")

# 3S (23 bp, 20 aligned): the read pointer must advance by the clip length, not by 1. The
# mismatch is at aligned offset 2 (bases of the clip do not count).
reader_json_test(softclip_shift softclip_shift.sam
  "summary.total_reads=1"
  "read1.length.reads.mapq_ge30.23=1"
  "read1.absolute.mismatches.mapq_ge30#len=3"
  "read1.absolute.mismatches.mapq_ge30.2=1")

# 5H10S40M and 40M10S5H: the length is the SEQ length (50, soft clip included), the position
# axes cover the 40 aligned bases
reader_json_test(hardclip_softclip hardclip_softclip.sam
  "summary.total_reads=2"
  "read1.length.reads.mapq_ge30.50=2"
  "read1.absolute.qual_bases.mapq_ge30#len=40")

# MD tags that disagree with the alignment: bounded by the alignment length
reader_json_test(md_inconsistent md_inconsistent.sam
  "summary.total_reads=3"
  "read1.length.reads.mapq_ge30.20=3"
  "read1.absolute.mismatches.mapq_ge30.10=1"
  "read1.absolute.mismatches.mapq_ge30#len=11")
reader_json_test(md_8digit md_8digit.sam
  "summary.total_reads=1"
  "read1.length.reads.mapq_ge30.20=1"
  "read1.absolute.mismatches.mapq_ge30#len=0")
reader_json_test(md_long_number md_long_number.sam
  "summary.total_reads=2"
  "read1.length.reads.mapq_ge30.20=2"
  "read1.absolute.mismatches.mapq_ge30.10=1")

# 300M reads: the 256 byte slack of the alignment buffers is too small to hide an
# unbounded MD walk (MD digits 99999 / 300A99999 against 300 aligned bases)
reader_json_test(md_long_alignment md_long_alignment.sam
  "summary.total_reads=3"
  "read1.length.reads.mapq_ge30.300=3"
  "read1.absolute.mismatches.mapq_ge30#len=0")

# NM and MD of the wrong type (written into the build tree): an MD that is not MD:Z is no MD
# (MD:i used to crash in strlen(NULL)), an NM that is not an integer is no NM (NM:Z:0 used to
# be read as 0, a perfect match). r1-r3: MD:i / MD:A / MD:B, no NM: no MD, identity NA.
# r4, r5: NM:Z / NM:f with MD:Z: MD counts, identity NA (10M has no = / X). r6: NM:i:1 with
# MD:Z:4C5, the control: identity 900.
set(READER_GEN_IN ${CMAKE_CURRENT_BINARY_DIR}/gen_in)
file(MAKE_DIRECTORY ${READER_GEN_IN})
set(aux_rec "\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\t")
file(WRITE ${READER_GEN_IN}/aux_types.sam
  "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
  "r1\t0${aux_rec}MD:i:4\n"
  "r2\t0${aux_rec}MD:A:x\n"
  "r3\t0${aux_rec}MD:B:c,1\n"
  "r4\t0${aux_rec}NM:Z:0\tMD:Z:10\n"
  "r5\t0${aux_rec}NM:f:0.5\tMD:Z:10\n"
  "r6\t0${aux_rec}NM:i:1\tMD:Z:4C5\n")
add_test(NAME reader:json:aux_types
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${READER_GEN_IN}/aux_types.sam
    -DOUTDIR=${READER_OUT}/aux_types
    "-DCHECKS=summary.total_reads=6\;summary.reads_without_md=3\;summary.identity_reads=1\;summary.identity_na=5"
    -P ${READER_HELPERS}/run_json.cmake)
set_tests_properties(reader:json:aux_types PROPERTIES LABELS reader)

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

reader_error_test(not_a_bam not_a_bam.bam "unrecognised input")
reader_error_test(truncated truncated.bam "truncated or corrupt")

# a 5000 bp read after three 20 bp reads: the length histogram grows with it (5000 bp is in the
# bin of edge 4953, bin 1157); the absolute position axis stops at 1000 bins and counts the other
# 4000 bases of the long read in bases_beyond_position_axis; the relative axis has 100 bins, the
# long read puts 50 bases into each (the short reads one base into every fifth bin).
reader_json_test(long_after_short long_after_short.sam
  "summary.total_reads=4"
  "summary.bases_beyond_position_axis=4000"
  "read1.length.reads.mapq_ge30.20=3"
  "read1.length.reads.mapq_ge30#len=1158"
  "read1.length.reads.mapq_ge30.1157=1"
  "read1.absolute.qual_bases.mapq_ge30.19=4"
  "read1.absolute.qual_bases.mapq_ge30.999=1"
  "read1.absolute.qual_bases.mapq_ge30#len=1000"
  "read1.relative.qual_bases.mapq_ge30#len=100"
  "read1.relative.qual_bases.mapq_ge30.0=53"
  "read1.relative.qual_bases.mapq_ge30.99=50")

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
  execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${dir} ${bam}
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

# --- WP3.2: -p cutoff counts -------------------------------------------------------------
# The reader stops at exactly N primary records, even in the middle of a batch, and counts
# secondary/supplementary records only up to that point. The input is written at test time
# (900 blocks of 10 records = 9000 records, more than two batches of 4096):
#   block = p p S p p U p p p p   (p primary, S secondary 0x100, U supplementary 0x800)
# With N = 8 b + r primaries (0 <= r < 8) the reader has seen b full blocks plus the records
# up to the r-th primary of the next block: secondary = b + (r >= 3), supplementary = b + (r >= 5).
# Each case runs with the default batch size and with 7 records per batch
# (SAMSTAT_BATCH_RECORDS, test-only): the counts must not depend on it.
# Cases are N:total_reads:secondary:supplementary:partial.
file(WRITE ${READER_GEN}/run_p_cutoff.cmake [=[
# -DSAMSTAT -DOUTDIR -DCHECK_JSON=<check_json.cmake> -DTHREADS=<n> -DCASES=<list>
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
set(block "")
foreach(rec "r1\t0" "r2\t0" "r3\t256" "r4\t0" "r5\t0" "r6\t2048" "r7\t0" "r8\t0" "r9\t0" "r10\t0")
  string(APPEND block "${rec}\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n")
endforeach()
string(REPEAT "${block}" 900 body)
set(input ${OUTDIR}/cutoff.sam)
file(WRITE ${input} "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n${body}")
if(NOT DEFINED THREADS)
  set(THREADS 1)
endif()
string(REPLACE "\;" ";" CASES "${CASES}")
foreach(batch default 7)
  if(batch STREQUAL "default")
    unset(ENV{SAMSTAT_BATCH_RECORDS})
  else()
    set(ENV{SAMSTAT_BATCH_RECORDS} ${batch})
  endif()
  foreach(case IN LISTS CASES)
    string(REPLACE ":" ";" f "${case}")
    list(GET f 0 n)
    list(GET f 1 total)
    list(GET f 2 sec)
    list(GET f 3 supp)
    list(GET f 4 partial)
    set(dir ${OUTDIR}/b${batch}_p${n})
    file(MAKE_DIRECTORY ${dir})
    execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -t ${THREADS} -p ${n} -d ${dir} ${input}
      RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
      message(FATAL_ERROR "-p ${n} (batch ${batch}): exit ${rc}\n${out}${err}")
    endif()
    execute_process(
      COMMAND ${CMAKE_COMMAND} -DJSON=${dir}/cutoff.sam.samstat.json
        "-DCHECKS=summary.total_reads=${total}\;summary.secondary=${sec}\;summary.supplementary=${supp}\;summary.is_partial_report=${partial}"
        -P ${CHECK_JSON}
      RESULT_VARIABLE rc2 OUTPUT_VARIABLE out2 ERROR_VARIABLE err2)
    if(NOT rc2 EQUAL 0)
      message(FATAL_ERROR "-p ${n} (batch ${batch}):\n${out2}${err2}")
    endif()
    file(REMOVE_RECURSE ${dir})
  endforeach()
endforeach()
]=])
# N=1, 6 (r=6, cutoff before the end of a block), 8 (a block boundary, one more primary
# follows), 4096 (a batch boundary), 4097, 5003 (N not a multiple of 4096), 7199 (the last
# primary is the next one), 7200 and 100000 (N is the whole file and more: not partial)
set(P_CUTOFF_CASES
  "1:1:0:0:true" "6:6:1:1:true" "8:8:1:1:true" "4096:4096:512:512:true"
  "4097:4097:512:512:true" "5003:5003:626:625:true" "7199:7199:900:900:true"
  "7200:7200:900:900:false" "100000:7200:900:900:false")
string(REPLACE ";" "\;" p_cutoff_cases "${P_CUTOFF_CASES}")
add_test(NAME reader:p_cutoff_counts
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${READER_OUT}/p_cutoff_counts
    -DCHECK_JSON=${READER_HELPERS}/check_json.cmake -DTHREADS=1 "-DCASES=${p_cutoff_cases}"
    -P ${READER_GEN}/run_p_cutoff.cmake)
set_tests_properties(reader:p_cutoff_counts PROPERTIES LABELS reader)

# --- CRAM (WP5.2) --------------------------------------------------------------------------
# The reference comes from REF_PATH / REF_CACHE only; the fixtures (tests/data/cram, see its
# README) carry M5 but no UR. REF_PATH=<dir>/%s finds the reference as <dir>/<md5>.
set(CRAM_DIR ${SAMSTAT_FIXTURE_DIR}/cram)
set(CRAM_REF_PATH "REF_PATH=${CRAM_DIR}/ref/%s")

# cram:with_ref  the same 20 reads as small.bam: the whole report equals the BAM golden (meta ignored)
add_test(NAME cram:with_ref
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${CRAM_DIR}/small.cram
    -DOUTDIR=${READER_OUT}/cram_with_ref
    -DGOLDEN=${PROJECT_SOURCE_DIR}/tests/golden/expected/small.bam.json
    "-DCHECKS=summary.total_reads=20\;summary.is_aligned=true"
    -P ${READER_HELPERS}/run_json.cmake)
set_tests_properties(cram:with_ref PROPERTIES LABELS reader ENVIRONMENT "${CRAM_REF_PATH}")

# cram:lzma_codec  a CRAM whose blocks are LZMA coded (needs the lzma library in htslib);
# 300 reads of 100 bp, 109 with MAPQ >= 30
add_test(NAME cram:lzma_codec
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${CRAM_DIR}/lzma.cram
    -DOUTDIR=${READER_OUT}/cram_lzma_codec
    "-DCHECKS=summary.total_reads=300\;summary.mapq_counts.read1.mapq_ge30=109\;summary.total_bases=30000\;read1.length.reads.mapq_ge30.100=109"
    -P ${READER_HELPERS}/run_json.cmake)
set_tests_properties(cram:lzma_codec PROPERTIES LABELS reader ENVIRONMENT "${CRAM_REF_PATH}")

# cram:no_ref  no reference reachable: exit 1, one message naming both variables, no report
# and no network access (libcurl is off in the bundled htslib); under a timeout.
# Registered for the bundled htslib only: a system htslib may be built with libcurl and
# then fetch the missing reference from the EBI server, so the test could reach the network.
if(NOT SAMSTAT_USE_SYSTEM_HTSLIB)
file(WRITE ${READER_GEN}/run_cram_no_ref.cmake [=[
# -DSAMSTAT -DINPUT=<cram> -DOUTDIR
cmake_minimum_required(VERSION 3.19)
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 1)
  message(FATAL_ERROR "expected exit 1, got ${rc}\nstdout:\n${out}\nstderr:\n${err}")
endif()
if(NOT "${out}${err}" MATCHES "REF_PATH[^\n]*REF_CACHE")
  message(FATAL_ERROR "the message does not name REF_PATH and REF_CACHE\nstdout:\n${out}\nstderr:\n${err}")
endif()
file(GLOB reports ${OUTDIR}/*)
if(reports)
  message(FATAL_ERROR "report written without a reference: ${reports}")
endif()
]=])
add_test(NAME cram:no_ref
  COMMAND ${CMAKE_COMMAND} -E env --unset=REF_PATH --unset=REF_CACHE
    ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${CRAM_DIR}/small.cram
    -DOUTDIR=${READER_OUT}/cram_no_ref -P ${READER_GEN}/run_cram_no_ref.cmake)
set_tests_properties(cram:no_ref PROPERTIES LABELS reader TIMEOUT 60)
endif()
