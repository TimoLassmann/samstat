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
  "summary.mapq_counts.read1.mapq_ge30=1"
  "summary.mapq_counts.read1.mapq_lt10=0"
  "read1.length.mapq_ge30.50=1")

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
