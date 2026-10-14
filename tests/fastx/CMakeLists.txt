# FASTA/FASTQ input tests (label: fastx).
#
# fastx:json:qual        FASTQ qualities are decoded with Phred offset 33 ('!' is Q0)
# fastx:json:<empty>     inputs without a sequence give an empty report, exit 0
# fastx:json:<no_nl>     a last record without trailing newline is not dropped
# fastx:warn:*           the empty-input warning is printed

set(FASTX_HELPERS ${PROJECT_SOURCE_DIR}/tests/cmake)
set(FASTX_OUT ${SAMSTAT_TEST_OUT_DIR}/fastx)

# fastx_json_test(<name> <fixture> <checks...>)
function(fastx_json_test name input)
  string(REPLACE ";" "\;" checks "${ARGN}")
  add_test(NAME fastx:json:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${SAMSTAT_FIXTURE_DIR}/${input}
      -DOUTDIR=${FASTX_OUT}/${name} "-DCHECKS=${checks}"
      -P ${FASTX_HELPERS}/run_json.cmake)
  set_tests_properties(fastx:json:${name} PROPERTIES LABELS fastx)
endfunction()

# qual.fq: r1 = ACGT / !III, r2 = ACGT / IIII ('!' = Q0, 'I' = Q40).
# With offset 0 the '!' would be counted in row 33 and the 'I' in row 41.
fastx_json_test(qual qual.fq
  "summary.total_reads=2"
  "summary.is_aligned=false"
  "read1.base_quality.unmapped.0.0=1"
  "read1.base_quality.unmapped.33.0=0"
  "read1.base_quality.unmapped.40.0=1"
  "read1.base_quality.unmapped.40.1=2"
  "read1.base_quality.unmapped.40.3=2"
  "read1.base_quality.unmapped.41.1=0")

# inputs that are valid but contain no sequence
fastx_json_test(empty_gz empty.fq.gz "summary.total_reads=0")
fastx_json_test(newline_only newline_only.fq "summary.total_reads=0")
fastx_json_test(header_only header_only.fa "summary.total_reads=0")

# the last record has no trailing newline
fastx_json_test(nonl_fq nonl.fq
  "summary.total_reads=1"
  "read1.length.unmapped.4=1"
  "read1.base_quality.unmapped.40.3=1")
fastx_json_test(nonl_fa nonl.fa
  "summary.total_reads=2"
  "read1.length.unmapped.4=2")

add_test(NAME fastx:warn:empty_gz
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DARGS=-d\;@OUTDIR@\;${SAMSTAT_FIXTURE_DIR}/empty.fq.gz"
    -DOUTDIR=${FASTX_OUT}/warn_empty_gz "-DEXPECT_STDERR_REGEX=no reads"
    "-DEXPECT_FILES=empty.fq.gz.samstat.html"
    -P ${FASTX_HELPERS}/run_expect.cmake)
set_tests_properties(fastx:warn:empty_gz PROPERTIES LABELS fastx)
