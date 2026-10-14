# JSON report tests (label: json).
#
# json:<file>      --json on every dev/*.sam|*.fa: parses, required keys, counts add up
# json:values:*    hand-verified values on small inputs (check_json.cmake)
# json:format_version_2  top-level layout of format_version 2
# json:no_v1_key:*  keys of format_version 1 and the "v2" wrapper are gone
# json:golden:*    full report compared with tests/golden/expected/<input>.json (meta ignored)
# json:off         without --json no .json file is written
# json:helper:*    self-tests of check_json.cmake / compare_golden.cmake
#
# Goldens live in tests/golden/expected/ and are regenerated in one way only:
#   cmake --build <build> --target update-goldens

set(JSON_HELPERS ${PROJECT_SOURCE_DIR}/tests/cmake)
set(JSON_OUT ${SAMSTAT_TEST_OUT_DIR}/json)
set(JSON_GOLDEN_DIR ${PROJECT_SOURCE_DIR}/tests/golden/expected)
set(JSON_HELPER_FIXTURE ${CMAKE_CURRENT_SOURCE_DIR}/data)

file(GLOB JSON_INPUTS ${SAMSTAT_TEST_DATA_DIR}/*.sam ${SAMSTAT_TEST_DATA_DIR}/*.fa)
foreach(input ${JSON_INPUTS})
  get_filename_component(name ${input} NAME)
  add_test(NAME json:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${input} -DOUTDIR=${JSON_OUT}/all_${name}
      -P ${JSON_HELPERS}/run_json.cmake)
  set_tests_properties(json:${name} PROPERTIES LABELS json)
endforeach()

# json_value_test(<name> <input> <checks...>)
function(json_value_test name input)
  # ';' inside an add_test argument has to be escaped to survive as one argument
  string(REPLACE ";" "\\;" checks "${ARGN}")
  add_test(NAME json:values:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${SAMSTAT_TEST_DATA_DIR}/${input}
      -DOUTDIR=${JSON_OUT}/values_${name} "-DCHECKS=${checks}"
      -P ${JSON_HELPERS}/run_json.cmake)
  set_tests_properties(json:values:${name} PROPERTIES LABELS json)
endfunction()

# json_golden_test(<name> <input>): compares with tests/golden/expected/<input>.json
function(json_golden_test name input)
  add_test(NAME json:golden:${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${SAMSTAT_TEST_DATA_DIR}/${input}
      -DOUTDIR=${JSON_OUT}/golden_${name} -DGOLDEN=${JSON_GOLDEN_DIR}/${input}.json
      -P ${JSON_HELPERS}/run_json.cmake)
  set_tests_properties(json:golden:${name} PROPERTIES LABELS json)
endfunction()

# xx#pair.sam: 3 read pairs, all MAPQ 1 (bin mapq_lt10), 10M each; read1 AAAAAAAAAA,
# read2 is stored reverse-complemented (TTTTTTTTTT -> AAAAAAAAAA); all qualities '*' = Q9
# (E[9] = round(10^-0.9 * 2^30) = 135176087, 3 bases per position = 405528261). No NM, MD or
# X operations: identity and mismatches are not available (reads_without_md = identity_na = 6).
# Relative axis: L = 10, bin = 10 * position, 3 bases in each of the bins 0, 10, ..., 90.
json_value_test(pair "xx#pair.sam"
  "format_version=2"
  "meta.input_file=xx#pair.sam"
  "summary.total_reads=6"
  "summary.read1_reads=3"
  "summary.read2_reads=3"
  "summary.paired=6"
  "summary.properly_paired=6"
  "summary.is_aligned=true"
  "summary.is_partial_report=false"
  "summary.total_bases=60"
  "summary.min_length=10"
  "summary.max_length=10"
  "summary.n50=10"
  "summary.median_length=10"
  "summary.mean_q>=8.999"
  "summary.mean_q<=9.001"
  "summary.pct_mapped>=99.99"
  "summary.pct_mapq_ge30<=0.001"
  "summary.median_identity="
  "summary.layout=short"
  "summary.reads_without_md=6"
  "summary.identity_reads=0"
  "summary.identity_na=6"
  "summary.mapq_counts.read1.mapq_lt10=3"
  "summary.mapq_counts.read1.mapq_ge30=0"
  "summary.mapq_counts.read2.mapq_lt10=3"
  "mapq_bins#len=4"
  "mapq_bins.0.id=mapq_ge30"
  "mapq_bins.3.id=unmapped"
  "mapq_bins.3.label=Unmapped"
  "axes.length_edges#len=1922"
  "axes.absolute.bins=1000"
  "axes.relative.bins=100"
  "axes.identity_bins=1001"
  "axes.quality_error_scale=1073741824"
  "read1.n_reads=3"
  "read1.length.reads.mapq_lt10#len=11"
  "read1.length.reads.mapq_lt10.10=3"
  "read1.length.bases.mapq_lt10.10=30"
  "read1.length.reads.mapq_ge30#len=0"
  "read1.read_mean_q.mapq_lt10#len=10"
  "read1.read_mean_q.mapq_lt10.9=3"
  "read1.identity.mapq_lt10#len=0"
  "read1.absolute.composition.mapq_lt10.A#len=10"
  "read1.absolute.composition.mapq_lt10.A.0=3"
  "read1.absolute.composition.mapq_lt10.A.9=3"
  "read1.absolute.composition.mapq_lt10.C#len=0"
  "read1.absolute.qual_bases.mapq_lt10#len=10"
  "read1.absolute.qual_bases.mapq_lt10.0=3"
  "read1.absolute.qual_error_sum.mapq_lt10.0=405528261"
  "read1.absolute.qual_error_sum.mapq_lt10.9=405528261"
  "read1.absolute.aligned.mapq_lt10.9=3"
  "read1.absolute.aligned_md.mapq_lt10#len=0"
  "read1.absolute.mismatches.mapq_lt10#len=0"
  "read1.absolute.insertions.mapq_lt10#len=0"
  "read1.relative.qual_bases.mapq_lt10#len=91"
  "read1.relative.qual_bases.mapq_lt10.0=3"
  "read1.relative.qual_bases.mapq_lt10.5=0"
  "read1.relative.qual_bases.mapq_lt10.90=3"
  "read2.n_reads=3"
  "read2.length.reads.mapq_lt10.10=3")

# xx.fa: two FASTA records of 20 nt (AAAAAAAAAATTTTTTTTTT), no alignment, no qualities.
json_value_test(fasta "xx.fa"
  "summary.total_reads=2"
  "summary.read1_reads=2"
  "summary.read2_reads=0"
  "summary.is_aligned=false"
  "summary.mean_q="
  "summary.pct_mapped="
  "summary.pct_mapq_ge30="
  "summary.pct_clipped="
  "summary.median_identity="
  "summary.mapq_counts.read1.unmapped=2"
  "read1.n_reads=2"
  "read1.length.reads.unmapped.20=2"
  "read1.length.bases.unmapped.20=40"
  "read1.absolute.composition.unmapped.A.0=2"
  "read1.absolute.composition.unmapped.A#len=10"
  "read1.absolute.composition.unmapped.T.10=2"
  "read1.absolute.composition.unmapped.T.19=2"
  "read1.absolute.qual_bases.unmapped#len=0"
  "read1.read_mean_q.unmapped#len=0")

# ce#unmap.sam: six unmapped reads (flag 4), no @SQ: is_aligned is false (has_sq or a mapped
# primary read; standard-plots.md 3.6).
json_value_test(unmapped "ce#unmap.sam"
  "summary.total_reads=6"
  "summary.is_aligned=false"
  "summary.mapq_counts.read1.unmapped=6"
  "summary.mapq_counts.read1.mapq_ge30=0"
  "summary.pct_mapped="
  "read1.n_reads=6"
  "read1.length.reads.unmapped.100=6"
  "read1.absolute.qual_bases.unmapped#len=100"
  "read1.absolute.qual_bases.unmapped.99=6")

# xx#blank.sam: header only, no reads.
json_value_test(blank "xx#blank.sam"
  "summary.total_reads=0"
  "summary.read1_reads=0"
  "read1.n_reads=0"
  "summary.min_length=0"
  "summary.n50=0"
  "read1.length.reads.unmapped#len=0")

json_golden_test(pair "xx#pair.sam")
json_golden_test(fasta "xx.fa")
json_golden_test(unmapped "ce#unmap.sam")
json_golden_test(ce2 "ce#2.sam")

add_test(NAME json:off
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${SAMSTAT_TEST_DATA_DIR}/xx#pair.sam"
    -DOUTDIR=${JSON_OUT}/off -DNO_JSON=ON
    -P ${JSON_HELPERS}/run_json.cmake)
set_tests_properties(json:off PROPERTIES LABELS json)

# --- helper self-tests -------------------------------------------------------

# check_json.cmake must fail on a wrong value, a missing path, a wrong length, a wrong boolean
set(JSON_BAD_CHECKS
  wrong_value "summary.total_reads=7" "mismatch for 'summary.total_reads': expected '7', got '6'"
  missing_path "summary.nope=1" "cannot resolve path 'summary.nope'"
  wrong_len "a_list#len=5" "mismatch for 'a_list': expected '5', got '3'"
  wrong_bool "flag=false" "mismatch for 'flag'")
list(LENGTH JSON_BAD_CHECKS n_bad)
math(EXPR last_bad "${n_bad}-1")
# PASS_REGULAR_EXPRESSION makes the test pass only when check_json fails for the
# stated reason (a missing mini.json or a syntax error in the script does not match).
foreach(i RANGE 0 ${last_bad} 3)
  math(EXPR j "${i}+1")
  math(EXPR k "${i}+2")
  list(GET JSON_BAD_CHECKS ${i} cname)
  list(GET JSON_BAD_CHECKS ${j} cchecks)
  list(GET JSON_BAD_CHECKS ${k} cregex)
  add_test(NAME json:helper:check_json_${cname}
    COMMAND ${CMAKE_COMMAND} -DJSON=${JSON_HELPER_FIXTURE}/mini.json "-DCHECKS=${cchecks}"
      -P ${JSON_HELPERS}/check_json.cmake)
  set_tests_properties(json:helper:check_json_${cname} PROPERTIES LABELS json PASS_REGULAR_EXPRESSION "${cregex}")
endforeach()

# an empty CHECKS list must fail instead of silently checking nothing
add_test(NAME json:helper:check_json_empty
  COMMAND ${CMAKE_COMMAND} -DJSON=${JSON_HELPER_FIXTURE}/mini.json "-DCHECKS=" -P ${JSON_HELPERS}/check_json.cmake)
set_tests_properties(json:helper:check_json_empty PROPERTIES LABELS json
  PASS_REGULAR_EXPRESSION "CHECKS contains no check")

string(REPLACE ";" "\\;" mini_ok "summary.total_reads=6;a_list#len=3;a_list.2=c;flag=true;summary.name=x y")
add_test(NAME json:helper:check_json_ok
  COMMAND ${CMAKE_COMMAND} -DJSON=${JSON_HELPER_FIXTURE}/mini.json "-DCHECKS=${mini_ok}"
    -P ${JSON_HELPERS}/check_json.cmake)
set_tests_properties(json:helper:check_json_ok PROPERTIES LABELS json)

# compare_golden.cmake ignores 'meta' but fails on any other difference
add_test(NAME json:helper:golden_ignores_meta
  COMMAND ${CMAKE_COMMAND} -DACTUAL=${JSON_HELPER_FIXTURE}/golden_a.json
    -DEXPECTED=${JSON_HELPER_FIXTURE}/golden_a_othermeta.json -P ${JSON_HELPERS}/compare_golden.cmake)
add_test(NAME json:helper:golden_detects_diff
  COMMAND ${CMAKE_COMMAND} -DACTUAL=${JSON_HELPER_FIXTURE}/golden_a.json
    -DEXPECTED=${JSON_HELPER_FIXTURE}/golden_b.json -P ${JSON_HELPERS}/compare_golden.cmake)
set_tests_properties(json:helper:golden_ignores_meta json:helper:golden_detects_diff PROPERTIES LABELS json)
# WILL_FAIL inverts the exit code; PASS_REGULAR_EXPRESSION then overrides it, so it
# checks the message names the first differing dotted path.
set_tests_properties(json:helper:golden_detects_diff PROPERTIES
  PASS_REGULAR_EXPRESSION "first differing path: summary.total_reads")

# --- schema ------------------------------------------------------------------
# json:schema_keys:<input>  member names of the report match docs/samstat.schema.json
foreach(input small.sam paired_small.sam empty.sam qual.fq supplementary.sam)
  add_test(NAME json:schema_keys:${input}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${SAMSTAT_FIXTURE_DIR}/${input}
      -DOUTDIR=${JSON_OUT}/schema_${input} -DSCHEMA=${PROJECT_SOURCE_DIR}/docs/samstat.schema.json
      -P ${JSON_HELPERS}/check_schema_keys.cmake)
  set_tests_properties(json:schema_keys:${input} PROPERTIES LABELS json)
endforeach()

# the version in meta is SAMStat's, not htslib's (config.h include order)
add_test(NAME json:meta_version
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DARGS=--json\;-d\;@OUTDIR@\;${SAMSTAT_FIXTURE_DIR}/small.sam"
    -DOUTDIR=${JSON_OUT}/meta_version -DJSON=small.sam.samstat.json
    "-DCHECKS=meta.samstat_version=${PROJECT_VERSION}"
    -P ${JSON_HELPERS}/run_expect.cmake)
set_tests_properties(json:meta_version PROPERTIES LABELS json)

# check_json.cmake >= / <= operators (mini.json: summary.total_reads=6, a_list has 3 elements)
string(REPLACE ";" "\\;" ops_ok "summary.total_reads>=6;summary.total_reads>=5;summary.total_reads<=6;summary.total_reads<=7;a_list#len>=3;a_list#len<=3")
add_test(NAME json:check_json_ops
  COMMAND ${CMAKE_COMMAND} -DJSON=${JSON_HELPER_FIXTURE}/mini.json "-DCHECKS=${ops_ok}"
    -P ${JSON_HELPERS}/check_json.cmake)
set_tests_properties(json:check_json_ops PROPERTIES LABELS json)

set(JSON_BAD_OPS
  ge_fails "summary.total_reads>=7" "mismatch for 'summary.total_reads': expected >= 7, got '6'"
  le_fails "summary.total_reads<=5" "mismatch for 'summary.total_reads': expected <= 5, got '6'"
  len_ge_fails "a_list#len>=4" "mismatch for 'a_list': expected >= 4, got '3'"
  non_numeric "summary.name>=1" "non-numeric operand")
list(LENGTH JSON_BAD_OPS n_bad_ops)
math(EXPR last_bad_ops "${n_bad_ops}-1")
foreach(i RANGE 0 ${last_bad_ops} 3)
  math(EXPR j "${i}+1")
  math(EXPR k "${i}+2")
  list(GET JSON_BAD_OPS ${i} cname)
  list(GET JSON_BAD_OPS ${j} cchecks)
  list(GET JSON_BAD_OPS ${k} cregex)
  add_test(NAME json:check_json_ops_${cname}
    COMMAND ${CMAKE_COMMAND} -DJSON=${JSON_HELPER_FIXTURE}/mini.json "-DCHECKS=${cchecks}"
      -P ${JSON_HELPERS}/check_json.cmake)
  set_tests_properties(json:check_json_ops_${cname} PROPERTIES LABELS json PASS_REGULAR_EXPRESSION "${cregex}")
endforeach()

# --- format_version 2 layout ---------------------------------------------------
# The v1 keys and the transitional "v2" wrapper are gone: the blocks of standard-plots.md
# section 7 are the top level. check_json.cmake only resolves paths, so absence is checked
# with PASS_REGULAR_EXPRESSION on "cannot resolve path".
add_test(NAME json:format_version_2
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${SAMSTAT_TEST_DATA_DIR}/xx#pair.sam"
    -DOUTDIR=${JSON_OUT}/format_version_2
    "-DCHECKS=format_version=2\\;summary.layout\\;axes.length_edges\\;read1.length.reads\\;read1.relative.composition.mapq_lt10.A"
    -P ${JSON_HELPERS}/run_json.cmake)
set_tests_properties(json:format_version_2 PROPERTIES LABELS json)
foreach(old read1.base_quality read1.base_composition read1.mismatches v2 read1.insertions summary.v2)
  add_test(NAME json:no_v1_key:${old}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> "-DINPUT=${SAMSTAT_TEST_DATA_DIR}/xx#pair.sam"
      -DOUTDIR=${JSON_OUT}/no_v1_${old} "-DCHECKS=${old}"
      -P ${JSON_HELPERS}/run_json.cmake)
  set_tests_properties(json:no_v1_key:${old} PROPERTIES LABELS json
    PASS_REGULAR_EXPRESSION "cannot resolve")
endforeach()
