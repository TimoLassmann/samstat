# MultiQC custom content output (label: multiqc).
#
#   multiqc:files        --multiqc writes the five <input>.samstat_<section>_mqc.json files
#                        (general, mapq, length, basequal, error; standard-plots.md section 8)
#   multiqc:off          without --multiqc no _mqc files are written
#   multiqc:<section>    structure and values of each file (ids, plot_type, hand-computed numbers)
#   multiqc:general_ids  exact list of general statistics columns for aligned, paired, FASTQ, FASTA
#   multiqc:error_absent_fastq  no error section for input that is not aligned
#   multiqc:relative_axis  both position sections have 100 x values for a short and a long sample
#   multiqc:run (optional) real MultiQC run on three samples, only if `multiqc` is on PATH

# check_json compares numbers as CMake normalises them: 50.00 reads as 50.0.
set(MQC_OUT ${SAMSTAT_TEST_OUT_DIR}/multiqc)
set(MQC_RUN ${PROJECT_SOURCE_DIR}/tests/cmake/run_mqc.cmake)
set(MQC_CHECK ${PROJECT_SOURCE_DIR}/tests/cmake/check_json.cmake)

# small.sam: 20 reads, 40M, 10x MAPQ 60 + 10x MAPQ 15, Q40 (tests/data/README.md)
# paired_small.sam: 2 pairs, R1 30M, R2 25M, MAPQ 60 and 20, Q40, flags 99/147

# One setup test produces all files so that the content checks can share them.
add_test(NAME multiqc:files
  COMMAND ${CMAKE_COMMAND}
    "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${SAMSTAT_FIXTURE_DIR}/small.sam"
    "-DOUTDIR=${MQC_OUT}/small" -DMULTIQC=ON -P ${MQC_RUN})
add_test(NAME multiqc:files_paired
  COMMAND ${CMAKE_COMMAND}
    "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${SAMSTAT_FIXTURE_DIR}/paired_small.sam"
    "-DOUTDIR=${MQC_OUT}/paired" -DMULTIQC=ON -P ${MQC_RUN})
add_test(NAME multiqc:off
  COMMAND ${CMAKE_COMMAND}
    "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${SAMSTAT_FIXTURE_DIR}/small.sam"
    "-DOUTDIR=${MQC_OUT}/off" -DMULTIQC=OFF -P ${MQC_RUN})
set_tests_properties(multiqc:files PROPERTIES FIXTURES_SETUP mqc_small)
set_tests_properties(multiqc:files_paired PROPERTIES FIXTURES_SETUP mqc_paired)

function(mqc_check name dir file)
  add_test(NAME multiqc:${name}
    COMMAND ${CMAKE_COMMAND} "-DJSON=${MQC_OUT}/${dir}/${file}" "-DCHECKS=${ARGN}" -P ${MQC_CHECK})
  if(dir STREQUAL "small")
    set_tests_properties(multiqc:${name} PROPERTIES FIXTURES_REQUIRED mqc_small)
  else()
    set_tests_properties(multiqc:${name} PROPERTIES FIXTURES_REQUIRED mqc_paired)
  endif()
  set_tests_properties(multiqc:${name} PROPERTIES LABELS multiqc)
endfunction()

# general stats: 20 reads, 10 with MAPQ >= 30 -> 50 %, none unmapped; unpaired -> no pair column;
# 40 bp reads at Q40 with MD: N50 40, mean Q 40, median identity 100. Percentages are checked with
# a tolerance because check_json compares numbers as CMake normalises them.
mqc_check(general small small.sam.samstat_general_mqc.json
  "id=samstat_general;plot_type=generalstats;pconfig#len=6;data.small.samstat_total_reads=20;data.small.samstat_pct_mapq_ge30>=49.99;data.small.samstat_pct_mapq_ge30<=50.01;data.small.samstat_pct_unmapped<=0.01;data.small.samstat_n50=40;data.small.samstat_mean_q>=39.99;data.small.samstat_mean_q<=40.01;data.small.samstat_median_identity>=99.9")
# paired: 4 reads, 2 MAPQ>=30, all properly paired
mqc_check(general_paired paired paired_small.sam.samstat_general_mqc.json
  "pconfig#len=7;data.paired_small.samstat_total_reads=4;data.paired_small.samstat_pct_mapq_ge30>=49.99;data.paired_small.samstat_pct_mapq_ge30<=50.01;data.paired_small.samstat_pct_proper_pair>=99.99;data.paired_small.samstat_n50=30")
mqc_check(mapq small small.sam.samstat_mapq_mqc.json
  "id=samstat_mapq;plot_type=table;parent_id=samstat;headers#len=4;data.small.mapq_ge30=10;data.small.mapq_10_29=10;data.small.mapq_lt10=0;data.small.unmapped=0")
# both read 1 and read 2 are summed in the table
mqc_check(mapq_paired paired paired_small.sam.samstat_mapq_mqc.json
  "data.paired_small.mapq_ge30=2;data.paired_small.mapq_10_29=2;data.paired_small.mapq_lt10=0")
mqc_check(length small small.sam.samstat_length_mqc.json
  "id=samstat_length;plot_type=linegraph;data.small#len=1;data.small.40=20")
# read length is read 1 only: 2 reads of length 30
mqc_check(length_paired paired paired_small.sam.samstat_length_mqc.json
  "data.paired_small.30=2;data.paired_small#len=1")
# Q40 everywhere: mean 40 in every populated bin of the relative axis (x = 1..100, the percent of
# the read length at which the bin ends). 40 bp reads fill 40 of the 100 bins, 30 bp reads 30.
mqc_check(basequal small small.sam.samstat_basequal_mqc.json
  "id=samstat_basequal;plot_type=linegraph;data.small#len=40;data.small.1>=39.99;data.small.1<=40.01;data.small.98>=39.99")
mqc_check(basequal_paired paired paired_small.sam.samstat_basequal_mqc.json
  "data.paired_small#len=30;data.paired_small.1>=39.99;data.paired_small.97>=39.99")
# error rate of the MAPQ >= 30 reads of read 1 (all bases match): 0 in every populated bin
mqc_check(error small small.sam.samstat_error_mqc.json
  "id=samstat_error;plot_type=linegraph;parent_id=samstat;data.small#len=40;data.small.1<=0.0001;data.small.98<=0.0001")
mqc_check(error_paired paired paired_small.sam.samstat_error_mqc.json
  "id=samstat_error;data.paired_small#len=30")

set_tests_properties(multiqc:files multiqc:files_paired multiqc:off PROPERTIES LABELS multiqc)

# Optional: a real MultiQC run. Registered only if multiqc is installed;
# -DSAMSTAT_REQUIRE_MULTIQC=ON turns a missing multiqc into a configure error (CI).
option(SAMSTAT_REQUIRE_MULTIQC "Fail the configure step if multiqc is not installed" OFF)
find_program(MULTIQC multiqc)
if(MULTIQC)
  set(MQC_INDIR ${MQC_OUT}/multiqc_in)
  add_test(NAME multiqc:run_setup_a
    COMMAND ${CMAKE_COMMAND} "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${SAMSTAT_FIXTURE_DIR}/small.sam"
      "-DOUTDIR=${MQC_INDIR}" -DMULTIQC=ON -P ${MQC_RUN})
  add_test(NAME multiqc:run_setup_b
    COMMAND ${CMAKE_COMMAND} "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${SAMSTAT_FIXTURE_DIR}/paired_small.sam"
      "-DOUTDIR=${MQC_INDIR}" -DMULTIQC=ON -DKEEP=ON -P ${MQC_RUN})
  add_test(NAME multiqc:run_setup_c
    COMMAND ${CMAKE_COMMAND} "-DSAMSTAT=$<TARGET_FILE:samstat>" "-DINPUT=${SAMSTAT_FIXTURE_DIR}/supplementary.sam"
      "-DOUTDIR=${MQC_INDIR}" -DMULTIQC=ON -DKEEP=ON -P ${MQC_RUN})
  add_test(NAME multiqc:run
    COMMAND ${CMAKE_COMMAND} "-DMULTIQC_EXE=${MULTIQC}" "-DINDIR=${MQC_INDIR}"
      "-DOUTDIR=${MQC_OUT}/multiqc_out" "-DSAMPLES=small;paired_small;supplementary"
      -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_multiqc_tool.cmake)
  # the three setups share one directory, so they must not run concurrently
  set_tests_properties(multiqc:run_setup_a multiqc:run_setup_b multiqc:run_setup_c PROPERTIES
    FIXTURES_SETUP mqc_tool RESOURCE_LOCK mqc_tool_dir LABELS multiqc)
  set_tests_properties(multiqc:run_setup_b PROPERTIES DEPENDS multiqc:run_setup_a)
  set_tests_properties(multiqc:run_setup_c PROPERTIES DEPENDS multiqc:run_setup_b)
  set_tests_properties(multiqc:run PROPERTIES FIXTURES_REQUIRED mqc_tool LABELS multiqc TIMEOUT 300)
  # the new error section is in the MultiQC report for every sample (all three are aligned, with MD)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/check_mqc_error_section.cmake [=[
# -DDATA=<multiqc_data dir> -DSAMPLES=<list>
set(f ${DATA}/multiqc_samstat_error_plot.txt)
if(NOT EXISTS ${f})
  message(FATAL_ERROR "MultiQC did not write ${f}")
endif()
file(READ ${f} content)
foreach(sample ${SAMPLES})
  string(FIND "${content}" "\n${sample}\t" pos)
  if(pos EQUAL -1)
    message(FATAL_ERROR "Sample '${sample}' missing in ${f}")
  endif()
endforeach()
]=])
  add_test(NAME multiqc:run_error_section
    COMMAND ${CMAKE_COMMAND} "-DDATA=${MQC_OUT}/multiqc_out/multiqc_data" "-DSAMPLES=small;paired_small;supplementary"
      -P ${CMAKE_CURRENT_BINARY_DIR}/check_mqc_error_section.cmake)
  set_tests_properties(multiqc:run_error_section PROPERTIES FIXTURES_REQUIRED mqc_tool DEPENDS multiqc:run LABELS multiqc)
  message(STATUS "multiqc found (${MULTIQC}): registering multiqc:run")
elseif(SAMSTAT_REQUIRE_MULTIQC)
  message(FATAL_ERROR "SAMSTAT_REQUIRE_MULTIQC is ON but multiqc was not found on PATH")
else()
  message(STATUS "multiqc not found: skipping optional multiqc:run test")
endif()

# --- escaping, sample names, multiple inputs, failure paths ------------------
# mqc_run_test(<name> <expect_rc> <args...>) - run_expect.cmake, label multiqc
function(mqc_run_test name rc)
  cmake_parse_arguments(T "" "REGEX;JSON;COPY_AS;FILE_TEXT;FILE_TEXT_IN" "ARGS;FILES;NO_FILES;MKDIR;CHECKS;COPY_FROM" ${ARGN})
  string(REPLACE ";" "\;" args "${T_ARGS}")
  string(REPLACE ";" "\;" files "${T_FILES}")
  string(REPLACE ";" "\;" checks "${T_CHECKS}")
  string(REPLACE ";" "\;" nofiles "${T_NO_FILES}")
  string(REPLACE ";" "\;" mkdir "${T_MKDIR}")
  set(cmd ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${MQC_OUT}/${name}
    "-DARGS=${args}" -DEXPECT_RC=${rc} "-DEXPECT_FILES=${files}" -DNO_NAN=ON)
  if(T_REGEX)
    list(APPEND cmd "-DEXPECT_STDERR_REGEX=${T_REGEX}")
  endif()
  if(T_NO_FILES)
    list(APPEND cmd "-DEXPECT_NO_FILES=${nofiles}")
  endif()
  if(T_MKDIR)
    list(APPEND cmd "-DMKDIR=${mkdir}")
  endif()
  if(T_COPY_FROM)
    list(APPEND cmd "-DCOPY_FROM=${T_COPY_FROM}" "-DCOPY_AS=${T_COPY_AS}")
  endif()
  if(T_JSON)
    list(APPEND cmd "-DJSON=${T_JSON}" "-DCHECKS=${checks}")
  endif()
  if(T_FILE_TEXT)
    list(APPEND cmd "-DFILE_TEXT=${T_FILE_TEXT}" "-DFILE_TEXT_IN=${T_FILE_TEXT_IN}")
  endif()
  list(APPEND cmd -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_expect.cmake)
  add_test(NAME multiqc:${name} COMMAND ${cmd})
  set_tests_properties(multiqc:${name} PROPERTIES LABELS multiqc)
endfunction()

# a file name with a double quote and a backslash must give valid JSON (check_json parses
# the file) whose sample key is the JSON-escaped name
mqc_run_test(escape_name 0
  COPY_FROM ${SAMSTAT_FIXTURE_DIR}/small.sam COPY_AS "a@Q@b@BS@c.sam"
  ARGS --multiqc -d @OUTDIR@ @COPY@
  FILES "a@Q@b@BS@c.sam.samstat_mapq_mqc.json"
  JSON "a@Q@b@BS@c.sam.samstat_mapq_mqc.json" CHECKS "id=samstat_mapq"
  FILE_TEXT "@Q@a@BS@@Q@b@BS@@BS@c@Q@: {"
  FILE_TEXT_IN "a@Q@b@BS@c.sam.samstat_mapq_mqc.json")

# sample name = file name without directory, '.gz' and one known extension
foreach(case "small.bam|small" "reads.fq.gz|reads" "READS.FQ|READS" "x.fa.gz|x" "y.fastq|y" "s1.CRAM|s1")
  string(REPLACE "|" ";" parts "${case}")
  list(GET parts 0 fname)
  list(GET parts 1 sample)
  if(fname MATCHES "^(small.bam|reads.fq.gz|READS.FQ)$")
    set(src ${SAMSTAT_FIXTURE_DIR}/${fname})
  elseif(fname STREQUAL "x.fa.gz")
    set(src ${SAMSTAT_FIXTURE_DIR}/reads.fa.gz)
  elseif(fname STREQUAL "y.fastq")
    set(src ${SAMSTAT_FIXTURE_DIR}/qual.fq)
  elseif(fname STREQUAL "s1.CRAM")
    # the type comes from the content (a SAM), the name is what is stripped: .cram in any case
    set(src ${SAMSTAT_FIXTURE_DIR}/small.sam)
  else()
    continue()
  endif()
  mqc_run_test(sample_${fname} 0
    COPY_FROM ${src} COPY_AS ${fname}
    ARGS --multiqc -d @OUTDIR@ @COPY@
    FILES ${fname}.samstat_mapq_mqc.json
    JSON ${fname}.samstat_mapq_mqc.json CHECKS "id=samstat_mapq" "data.${sample}")
endforeach()

# several inputs and threads: four files per input, one distinct sample key each
foreach(threads 1 4)
  add_test(NAME multiqc:multi_input_t${threads}
    COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${MQC_OUT}/multi_t${threads}
      "-DARGS=-t\;${threads}\;--multiqc\;-d\;@OUTDIR@\;${SAMSTAT_FIXTURE_DIR}/small.sam\;${SAMSTAT_FIXTURE_DIR}/paired_small.sam\;${SAMSTAT_FIXTURE_DIR}/small.bam"
      "-DEXPECT_FILES=small.sam.samstat_mapq_mqc.json\;paired_small.sam.samstat_mapq_mqc.json\;small.bam.samstat_mapq_mqc.json\;small.sam.samstat_general_mqc.json\;paired_small.sam.samstat_length_mqc.json\;small.bam.samstat_basequal_mqc.json"
      -DJSON=paired_small.sam.samstat_mapq_mqc.json "-DCHECKS=data.paired_small.mapq_ge30=2"
      -P ${PROJECT_SOURCE_DIR}/tests/cmake/run_expect.cmake)
  set_tests_properties(multiqc:multi_input_t${threads} PROPERTIES LABELS multiqc)
endforeach()

# an input without reads: general and mapq files (valid, no nan/inf); length and
# base quality are not written (docs/multiqc.md)
mqc_run_test(empty_input 0
  ARGS --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/empty.sam
  FILES empty.sam.samstat_general_mqc.json empty.sam.samstat_mapq_mqc.json
  NO_FILES empty.sam.samstat_length_mqc.json empty.sam.samstat_basequal_mqc.json
  JSON empty.sam.samstat_general_mqc.json CHECKS "id=samstat_general" "data.empty.samstat_total_reads=0")

# --json and --multiqc together write both
mqc_run_test(with_json 0
  ARGS --json --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/small.sam
  FILES small.sam.samstat.json small.sam.samstat_general_mqc.json small.sam.samstat.html
  JSON small.sam.samstat.json CHECKS "summary.total_reads=20")

# an output directory that does not exist is created (mkdir -p)
mqc_run_test(missing_dir 0
  ARGS --multiqc -d @OUTDIR@/does/not/exist ${SAMSTAT_FIXTURE_DIR}/small.sam
  FILES does/not/exist/small.sam.samstat_mapq_mqc.json)

# a MultiQC file that cannot be created (a directory has its name): exit 1 with a message
mqc_run_test(unwritable_file 1
  MKDIR small.sam.samstat_mapq_mqc.json
  ARGS --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/small.sam
  REGEX "Failed to open output file")

# --- standard-plots.md section 8 ----------------------------------------------------------
# mqc_sample_check(<name> <input> <file name = sample.ext> <section> <checks...>): run samstat
# --multiqc on <input> in its own directory and check one section file
function(mqc_sample_check name input fname section)
  mqc_run_test(${name} 0
    ARGS --multiqc -d @OUTDIR@ ${input}
    FILES ${fname}.samstat_${section}_mqc.json
    JSON ${fname}.samstat_${section}_mqc.json CHECKS ${ARGN})
endfunction()

# General statistics: the exact id list per kind of input (pconfig is an array of one-key objects
# and data.<sample> holds the same ids). A column that does not apply is absent.
# aligned and unpaired
mqc_check(general_ids small small.sam.samstat_general_mqc.json
  "pconfig#len=6;pconfig.0.samstat_total_reads;pconfig.1.samstat_pct_mapq_ge30;pconfig.2.samstat_pct_unmapped;pconfig.3.samstat_n50;pconfig.4.samstat_mean_q;pconfig.5.samstat_median_identity;data.small#len=6")
# aligned and paired: the properly paired column follows the unmapped column
mqc_check(general_ids_paired paired paired_small.sam.samstat_general_mqc.json
  "pconfig#len=7;pconfig.0.samstat_total_reads;pconfig.1.samstat_pct_mapq_ge30;pconfig.2.samstat_pct_unmapped;pconfig.3.samstat_pct_proper_pair;pconfig.4.samstat_n50;pconfig.5.samstat_mean_q;pconfig.6.samstat_median_identity;data.paired_small#len=7")
# FASTQ (qual.fq: 2 reads of 4 bp with qualities): no alignment columns
mqc_sample_check(general_ids_fastq ${SAMSTAT_FIXTURE_DIR}/qual.fq qual.fq general
  "pconfig#len=3;pconfig.0.samstat_total_reads;pconfig.1.samstat_n50;pconfig.2.samstat_mean_q;data.qual#len=3;data.qual.samstat_total_reads=2;data.qual.samstat_n50=4")
# FASTA (reads.fa.gz: 20, 25 and 30 bp, N50 25): no alignment columns and no mean Q
mqc_sample_check(general_ids_fasta ${SAMSTAT_FIXTURE_DIR}/reads.fa.gz reads.fa.gz general
  "pconfig#len=2;pconfig.0.samstat_total_reads;pconfig.1.samstat_n50;data.reads#len=2;data.reads.samstat_total_reads=3;data.reads.samstat_n50=25")
# no reads: only the read count
mqc_sample_check(general_ids_empty ${SAMSTAT_FIXTURE_DIR}/empty.sam empty.sam general
  "pconfig#len=1;data.empty#len=1;data.empty.samstat_total_reads=0")

# No error section for input that is not aligned: FASTQ, FASTA and a SAM without @SQ and
# without mapped reads. The other sections are written.
mqc_run_test(error_absent_fastq 0
  ARGS --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/qual.fq
  FILES qual.fq.samstat_general_mqc.json qual.fq.samstat_basequal_mqc.json
  NO_FILES qual.fq.samstat_error_mqc.json)
mqc_run_test(error_absent_fasta 0
  ARGS --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/reads.fa.gz
  FILES reads.fa.gz.samstat_general_mqc.json reads.fa.gz.samstat_length_mqc.json
  NO_FILES reads.fa.gz.samstat_error_mqc.json reads.fa.gz.samstat_basequal_mqc.json)
mqc_run_test(error_absent_ubam 0
  ARGS --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/longread/ubam_like.sam
  FILES ubam_like.sam.samstat_general_mqc.json ubam_like.sam.samstat_basequal_mqc.json
  NO_FILES ubam_like.sam.samstat_error_mqc.json)

# Hand-computed N50 and median identity (fixtures of tests/data/longread, see the README there).
# n50_five.sam: reads of 100, 200, 300, 400 and 600 bp, Q30, identity 100 %: N50 400
mqc_sample_check(n50_five ${SAMSTAT_FIXTURE_DIR}/longread/n50_five.sam n50_five.sam general
  "data.n50_five.samstat_n50=400;data.n50_five.samstat_median_identity>=99.9;data.n50_five.samstat_mean_q>=29.99;data.n50_five.samstat_mean_q<=30.01")
# identity_indel.sam: identity bins 980, 989, 990 and 1000 (tenths of a percent); the median is the
# lower edge of the second bin, 98.9 %
mqc_sample_check(median_identity ${SAMSTAT_FIXTURE_DIR}/longread/identity_indel.sam identity_indel.sam general
  "data.identity_indel.samstat_median_identity>=98.85;data.identity_indel.samstat_median_identity<=98.95;data.identity_indel.samstat_n50=100")
# display length bins: 1 bp exact up to 100, then edges 101, 103, ... (about 1 % wide): the reads
# of 100, 200, 300, 400 and 600 bp are in the bins with lower edges 100, 199, 300, 399 and 599
mqc_sample_check(length_display_bins ${SAMSTAT_FIXTURE_DIR}/longread/n50_five.sam n50_five.sam length
  "pconfig.xLog=true;data.n50_five#len=5;data.n50_five.100=1;data.n50_five.199=1;data.n50_five.300=1;data.n50_five.399=1;data.n50_five.599=1")

# Both position sections always use the relative axis: 100 x values for a short and a long
# sample (3 reads of 150 bp and of 1500 bp, all bases match, Q40). The inputs are written into
# the build tree at configure time.
set(MQC_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${MQC_GEN})
foreach(len 150 1500)
  math(EXPR quarter "${len}/4+1")
  string(REPEAT "ACGT" ${quarter} seq)
  string(SUBSTRING "${seq}" 0 ${len} seq)
  string(REPEAT "I" ${len} qual)
  set(sam "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n")
  foreach(i 1 2 3)
    string(APPEND sam "r${i}\t0\tchr1\t${i}\t60\t${len}M\t*\t0\t0\t${seq}\t${qual}\tNM:i:0\tMD:Z:${len}\n")
  endforeach()
  file(WRITE ${MQC_GEN}/len${len}.sam "${sam}")
  foreach(section basequal error)
    mqc_run_test(relative_axis:${len}:${section} 0
      ARGS --multiqc -d @OUTDIR@ ${MQC_GEN}/len${len}.sam
      FILES len${len}.sam.samstat_${section}_mqc.json
      JSON len${len}.sam.samstat_${section}_mqc.json
      CHECKS "data.len${len}#len=100" "data.len${len}.1" "data.len${len}.100"
             "pconfig.xlab=position in read (% of read length)" "pconfig.xmin=1" "pconfig.xmax=100")
  endforeach()
endforeach()
