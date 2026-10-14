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

