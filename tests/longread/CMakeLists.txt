# Long-read and summary value tests (label: longread).
#
# The expected values come from tests/data/longread/README.md, i.e. from the fixture generator's ground truth
# and the definitions of docs/plan-v2.4/standard-plots.md section 3; none was read from a
# samstat run. Keys are the `summary` keys of standard-plots.md section 7.
#
# Transitional prefix: in waves 3-4 the v2 blocks live under a top-level "v2" key next to
# the v1 keys. WP4.3 promotes them to the top level and sets LR_PREFIX to the empty string
# (and passes -DPREFIX= to run_ultralong.cmake) in this one place.
#
# Tolerances (testing T3): counts exact; N50 and median length within 1%; median identity
# within 0.5 percentage points; mean Q within 0.1; other percentages within 0.01.
#
# longread:value:<fixture>   summary values of one fixture
# longread:mapq_sum:<fixture> MAPQ class counts add up to the primary reads
# longread:ultralong          a 1 Mb read: JSON < 1 MB and bounded array lengths

set(LR_PREFIX "v2")
set(LR_DIR ${SAMSTAT_FIXTURE_DIR}/longread)
set(LR_OUT ${SAMSTAT_TEST_OUT_DIR}/longread)
set(LR_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${LR_GEN})

if(LR_PREFIX STREQUAL "")
  set(LR_P "")
else()
  set(LR_P "${LR_PREFIX}.")
endif()

# --- decimal helpers (CMake has integer arithmetic only) -------------------------------------
# lr_dec2int(<out> <decimal string>): value x 10000 ("92.9069" -> 929069, "100" -> 1000000)
function(lr_dec2int out s)
  if(NOT s MATCHES "^([0-9]+)(\\.([0-9]+))?$")
    message(FATAL_ERROR "longread tests: bad decimal '${s}'")
  endif()
  set(ip "${CMAKE_MATCH_1}")
  set(fp "${CMAKE_MATCH_3}0000")
  string(SUBSTRING "${fp}" 0 4 fp)
  string(REGEX REPLACE "^0+([0-9])" "\\1" fp "${fp}")
  string(REGEX REPLACE "^0+([0-9])" "\\1" ip "${ip}")
  math(EXPR v "${ip}*10000+${fp}")
  set(${out} ${v} PARENT_SCOPE)
endfunction()
# lr_int2dec(<out> <value x 10000>): 929069 -> "92.9069"
function(lr_int2dec out v)
  if(v LESS 0)
    set(v 0)
  endif()
  math(EXPR ip "${v}/10000")
  math(EXPR fp "${v}%10000")
  string(LENGTH "${fp}" n)
  math(EXPR pad "4-${n}")
  if(pad GREATER 0)
    string(REPEAT "0" ${pad} zeros)
    set(fp "${zeros}${fp}")
  endif()
  set(${out} "${ip}.${fp}" PARENT_SCOPE)
endfunction()

# lr_value_test(<fixture> [EQ k=v...] [PCT1 k=exact...] [NEAR k=center:tolerance...])
#   EQ    exact: integers, booleans (true / false), strings (long / short), null
#   PCT1  integer within 1% of the exact value (>= 0.99 x and <= 1.01 x, rounded outwards)
#   NEAR  decimal within +-tolerance of the center
function(lr_value_test fixture)
  cmake_parse_arguments(T "" "" "EQ;PCT1;NEAR" ${ARGN})
  set(checks "")
  foreach(kv IN LISTS T_EQ)
    # check_json reads a JSON null as the empty string (string(JSON GET)): "k=null" -> "k="
    string(REGEX REPLACE "=null$" "=" kv "${kv}")
    list(APPEND checks "${LR_P}summary.${kv}")
  endforeach()
  foreach(kv IN LISTS T_PCT1)
    string(REPLACE "=" ";" parts "${kv}")
    list(GET parts 0 key)
    list(GET parts 1 exact)
    math(EXPR lo "${exact}*99/100")
    math(EXPR hi "(${exact}*101+99)/100")
    list(APPEND checks "${LR_P}summary.${key}>=${lo}" "${LR_P}summary.${key}<=${hi}")
  endforeach()
  foreach(kv IN LISTS T_NEAR)
    string(REGEX MATCH "^([a-z0-9_]+)=([0-9.]+):([0-9.]+)$" m "${kv}")
    if(NOT m)
      message(FATAL_ERROR "longread tests: bad NEAR spec '${kv}'")
    endif()
    set(key ${CMAKE_MATCH_1})
    lr_dec2int(c "${CMAKE_MATCH_2}")
    lr_dec2int(tol "${CMAKE_MATCH_3}")
    math(EXPR lo "${c}-${tol}")
    math(EXPR hi "${c}+${tol}")
    lr_int2dec(lo_s ${lo})
    lr_int2dec(hi_s ${hi})
    list(APPEND checks "${LR_P}summary.${key}>=${lo_s}" "${LR_P}summary.${key}<=${hi_s}")
  endforeach()
  string(REPLACE ";" "\;" checks_arg "${checks}")
  add_test(NAME longread:value:${fixture}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${LR_DIR}/${fixture}
      -DOUTDIR=${LR_OUT}/value_${fixture} "-DCHECKS=${checks_arg}"
      -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_json.cmake)
  set_tests_properties(longread:value:${fixture} PROPERTIES LABELS longread)
endfunction()

# Single-ended fixtures put all reads in read1: the MAPQ class counts are checked as
# summary.mapq_counts.read1.<id>. `mq` builds the four EQ items.
function(lr_mq out ge30 mid low unmapped)
  set(${out}
    "mapq_counts.read1.mapq_ge30=${ge30}" "mapq_counts.read1.mapq_10_29=${mid}"
    "mapq_counts.read1.mapq_lt10=${low}" "mapq_counts.read1.unmapped=${unmapped}"
    PARENT_SCOPE)
endfunction()

# --- ont_like.sam: 24 reads, 22 mapped, long layout ------------------------------------------
lr_mq(mq 12 4 6 2)
lr_value_test(ont_like.sam
  EQ total_reads=24 read1_reads=24 read2_reads=0 paired=0 properly_paired=0 secondary=0
     supplementary=0 is_aligned=true total_bases=36990 min_length=538 max_length=3442
     layout=long identity_reads=22 identity_na=0 reads_without_md=0 ${mq}
  PCT1 n50=2090 median_length=1283
  NEAR mean_q=10.4577:0.1 pct_mapped=91.6667:0.01 pct_mapq_ge30=50:0.01 pct_clipped=5.5942:0.01
       median_identity=92.9069:0.5)

# --- ont_like.bam: the same plus a 52.5 kb read with a 70000-operation CIGAR (CG tag) ----------
lr_mq(mq 13 4 6 2)
lr_value_test(ont_like.bam
  EQ total_reads=25 read1_reads=25 read2_reads=0 is_aligned=true total_bases=89490
     min_length=538 max_length=52500 layout=long identity_reads=23 identity_na=0
     reads_without_md=0 ${mq}
  PCT1 n50=52500 median_length=1434
  NEAR mean_q=13.6586:0.1 pct_mapped=92:0.01 pct_mapq_ge30=52:0.01 pct_clipped=2.167:0.01
       median_identity=92.9069:0.5)

# --- supp_chain.sam: secondary / supplementary records are not reads -------------------------
lr_mq(mq 2 0 0 0)
lr_value_test(supp_chain.sam
  EQ total_reads=2 read1_reads=2 secondary=1 supplementary=2 is_aligned=true total_bases=3000
     min_length=1000 max_length=2000 layout=long identity_reads=2 identity_na=0 ${mq}
  PCT1 n50=2000 median_length=1000
  NEAR mean_q=15.1028:0.1 pct_mapped=100:0.01 pct_mapq_ge30=100:0.01
       pct_clipped=33.3333:0.01 median_identity=99.9:0.5)

# --- hifi_like.sam: high qualities do not saturate ---------------------------------------------
lr_mq(mq 3 1 0 0)
lr_value_test(hifi_like.sam
  EQ total_reads=4 read1_reads=4 is_aligned=true total_bases=18000 min_length=3000
     max_length=6000 layout=long identity_reads=4 identity_na=0 ${mq}
  PCT1 n50=4800 median_length=4200
  NEAR mean_q=37.0198:0.1 pct_mapped=100:0.01 pct_mapq_ge30=75:0.01 pct_clipped=0:0.01
       median_identity=99.7001:0.5)

# --- long.fq: mean Q over error probabilities (11.89), not over Phred values (14.98) --------------
lr_value_test(long.fq
  EQ total_reads=12 read1_reads=12 is_aligned=false total_bases=19391 min_length=1
     max_length=5000 layout=long mapq_counts.read1.unmapped=12
  PCT1 n50=2600 median_length=1200
  NEAR mean_q=11.8917:0.1)

# --- eqx_nomd.sam: =/X operations without MD; reads 4-5 without NM use the X count ------------
lr_mq(mq 4 1 0 0)
lr_value_test(eqx_nomd.sam
  EQ total_reads=5 read1_reads=5 is_aligned=true total_bases=2530 min_length=300
     max_length=700 layout=short identity_reads=5 identity_na=0 reads_without_md=0 ${mq}
  PCT1 n50=600 median_length=510
  NEAR mean_q=13.8782:0.1 pct_mapq_ge30=80:0.01 pct_clipped=2.6087:0.01
       median_identity=95.082:0.5)

# --- nonm.sam: MD but no NM and plain M: identity not available ----------------------------------
lr_mq(mq 4 0 1 0)
lr_value_test(nonm.sam
  EQ total_reads=5 read1_reads=5 is_aligned=true total_bases=2300 layout=short
     identity_reads=0 identity_na=5 reads_without_md=0 median_identity=null ${mq}
  PCT1 n50=520 median_length=450
  NEAR mean_q=15.0286:0.1 pct_mapq_ge30=80:0.01 pct_clipped=0:0.01)

# --- ubam_like.sam: no @SQ, every record unmapped: unaligned layout -------------------------------
lr_value_test(ubam_like.sam
  EQ total_reads=6 total_bases=4100 min_length=300 max_length=1100
     layout=short mapq_counts.read1.unmapped=6
  PCT1 n50=750 median_length=600
  NEAR mean_q=12.1286:0.1)

# is_aligned=false for a uBAM without @SQ needs the reader to supply has_sq (WP3.2 removes
# the input-type assignment in samstat.c). Until WP3.2 is merged (gate G3) this check stays
# disabled; the collector side is covered by collect:unit.
add_test(NAME longread:is_aligned:ubam_like.sam
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${LR_DIR}/ubam_like.sam
    -DOUTDIR=${LR_OUT}/is_aligned_ubam_like "-DCHECKS=${LR_P}summary.is_aligned=false"
    -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_json.cmake)
set_tests_properties(longread:is_aligned:ubam_like.sam PROPERTIES LABELS longread DISABLED TRUE)   # enable at G3, after WP3.2

# --- sq_nomap.sam: @SQ present, no mapped read: aligned layout, 100% unmapped ---------------------
lr_mq(mq 0 0 0 6)
lr_value_test(sq_nomap.sam
  EQ total_reads=6 read1_reads=6 is_aligned=true total_bases=3820 layout=short
     identity_reads=0 identity_na=0 pct_clipped=null median_identity=null ${mq}
  PCT1 n50=700 median_length=550
  NEAR mean_q=12.0356:0.1 pct_mapped=0:0.01 pct_mapq_ge30=0:0.01)

# --- mapq255.sam: MAPQ 255 counts as >= 30 ------------------------------------------------------
lr_mq(mq 6 0 1 0)
lr_value_test(mapq255.sam
  EQ total_reads=7 read1_reads=7 is_aligned=true total_bases=3850 min_length=320
     max_length=760 layout=short identity_reads=7 identity_na=0 ${mq}
  PCT1 n50=640 median_length=560
  NEAR mean_q=28.0353:0.1 pct_mapped=100:0.01 pct_mapq_ge30=85.7143:0.01 pct_clipped=0:0.01
       median_identity=96.7235:0.5)

# --- single_len.sam: 10 proper pairs of 150 bp reads (short layout, axis ends at 150) --------------
lr_value_test(single_len.sam
  EQ total_reads=20 read1_reads=10 read2_reads=10 paired=20 properly_paired=20
     is_aligned=true total_bases=3000 min_length=150 max_length=150 layout=short
     identity_reads=20 identity_na=0
  PCT1 n50=150 median_length=150
  NEAR mean_q=30.0519:0.1 pct_mapped=100:0.01 pct_mapq_ge30=70:0.01 pct_clipped=0:0.01
       median_identity=99.3333:0.5)

# --- hand-computable cases ---------------------------------------------------------------------
# N50 of 100, 200, 300, 400, 600 = 400; median 300
lr_value_test(n50_five.sam
  EQ total_reads=5 total_bases=1600 layout=short
  PCT1 n50=400 median_length=300
  NEAR median_identity=100:0.5)
# gap-compressed identity (M - X) / (M + I events + D events): bins 980, 989, 990, 1000
lr_value_test(identity_indel.sam
  EQ total_reads=4 total_bases=400 identity_reads=4 identity_na=0
  NEAR median_identity=98.9796:0.5)
# soft-clipped bases of mapped primaries: 60 of 400 = 15%; one unmapped read stays in total_bases
lr_value_test(softclip_pct.sam
  EQ total_reads=5 total_bases=500 mapq_counts.read1.mapq_ge30=4 mapq_counts.read1.unmapped=1
  NEAR pct_mapped=80:0.01 pct_mapq_ge30=80:0.01 pct_clipped=15:0.01)

# --- MAPQ classes add up to the primary reads ------------------------------------------------------
# -DEXPECT=<ge30>,<10_29>,<lt10>,<unmapped>: class totals over both mates
file(WRITE ${LR_GEN}/check_mapq_sum.cmake [=[
# -DSAMSTAT -DINPUT -DOUTDIR -DPREFIX=<v2 object path, may be empty> -DEXPECT=a,b,c,d
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
execute_process(COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} on ${INPUT}\n${out}\n${err}")
endif()
get_filename_component(name ${INPUT} NAME)
file(READ ${OUTDIR}/${name}.samstat.json doc)
if(PREFIX STREQUAL "")
  set(base "")
else()
  set(base ${PREFIX})
endif()
string(JSON total GET "${doc}" ${base} summary total_reads)
string(JSON r2 GET "${doc}" ${base} summary read2_reads)
set(mates read1)
if(r2 GREATER 0)
  list(APPEND mates read2)
endif()
set(classes mapq_ge30 mapq_10_29 mapq_lt10 unmapped)
string(REPLACE "," ";" expect "${EXPECT}")
set(sum 0)
set(failed "")
set(i 0)
foreach(c IN LISTS classes)
  set(n 0)
  foreach(m IN LISTS mates)
    string(JSON v GET "${doc}" ${base} summary mapq_counts ${m} ${c})
    math(EXPR n "${n}+${v}")
  endforeach()
  list(GET expect ${i} want)
  if(NOT n EQUAL want)
    string(APPEND failed "  ${c}: ${n}, expected ${want}\n")
  endif()
  math(EXPR sum "${sum}+${n}")
  math(EXPR i "${i}+1")
endforeach()
if(NOT sum EQUAL total)
  string(APPEND failed "  classes add up to ${sum}, total_reads is ${total}\n")
endif()
if(failed)
  message(FATAL_ERROR "${name}: MAPQ classes\n${failed}")
endif()
]=])
function(lr_mapq_sum_test fixture expect)
  add_test(NAME longread:mapq_sum:${fixture}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${LR_DIR}/${fixture}
      -DOUTDIR=${LR_OUT}/mapq_sum_${fixture} -DPREFIX=${LR_PREFIX} -DEXPECT=${expect}
      -P ${LR_GEN}/check_mapq_sum.cmake)
  set_tests_properties(longread:mapq_sum:${fixture} PROPERTIES LABELS longread)
endfunction()
lr_mapq_sum_test(ont_like.sam 12,4,6,2)
lr_mapq_sum_test(ont_like.bam 13,4,6,2)
lr_mapq_sum_test(single_len.sam 14,4,2,0)
lr_mapq_sum_test(supp_chain.sam 2,0,0,0)
lr_mapq_sum_test(mapq255.sam 6,0,1,0)
lr_mapq_sum_test(sq_nomap.sam 0,0,0,6)

# --- T4: a 1 Mb read gives a small JSON with bounded arrays ------------------------------------------
add_test(NAME longread:ultralong
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${LR_OUT}/ultralong -DPREFIX=${LR_PREFIX}
    -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_ultralong.cmake)
set_tests_properties(longread:ultralong PROPERTIES LABELS longread)
