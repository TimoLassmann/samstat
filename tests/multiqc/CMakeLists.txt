# MultiQC custom content output (label: multiqc).
#
#   multiqc:files        --multiqc writes the four <input>.samstat_<section>_mqc.json files
#   multiqc:off          without --multiqc no _mqc files are written
#   multiqc:<section>    structure and values of each file (ids, plot_type, hand-computed numbers)
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

# general stats: 20 reads, 10 with MAPQ >= 30 -> 50 %, none unmapped; unpaired -> no pair column
mqc_check(general small small.sam.samstat_general_mqc.json
  "id=samstat_general;plot_type=generalstats;pconfig#len=3;data.small.samstat_total_reads=20;data.small.samstat_pct_mapq_ge30=50.0;data.small.samstat_pct_unmapped=0.0")
# paired: 4 reads, 2 MAPQ>=30, all properly paired
mqc_check(general_paired paired paired_small.sam.samstat_general_mqc.json
  "pconfig#len=4;data.paired_small.samstat_total_reads=4;data.paired_small.samstat_pct_mapq_ge30=50.0;data.paired_small.samstat_pct_proper_pair=100.0")
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
# Q40 everywhere: mean 40 at all positions
mqc_check(basequal small small.sam.samstat_basequal_mqc.json
  "id=samstat_basequal;plot_type=linegraph;data.small#len=40;data.small.1=40.0;data.small.40=40.0")
mqc_check(basequal_paired paired paired_small.sam.samstat_basequal_mqc.json
  "data.paired_small#len=30;data.paired_small.30=40.0")

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
foreach(case "small.bam|small" "reads.fq.gz|reads" "READS.FQ|READS" "x.fa.gz|x" "y.fastq|y")
  string(REPLACE "|" ";" parts "${case}")
  list(GET parts 0 fname)
  list(GET parts 1 sample)
  if(fname MATCHES "^(small.bam|reads.fq.gz|READS.FQ)$")
    set(src ${SAMSTAT_FIXTURE_DIR}/${fname})
  elseif(fname STREQUAL "x.fa.gz")
    set(src ${SAMSTAT_FIXTURE_DIR}/reads.fa.gz)
  elseif(fname STREQUAL "y.fastq")
    set(src ${SAMSTAT_FIXTURE_DIR}/qual.fq)
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

# an output directory that does not exist is a usage error
mqc_run_test(missing_dir 2
  ARGS --multiqc -d @OUTDIR@/does/not/exist ${SAMSTAT_FIXTURE_DIR}/small.sam
  REGEX "does not exist")

# a MultiQC file that cannot be created (a directory has its name): exit 1 with a message
mqc_run_test(unwritable_file 1
  MKDIR small.sam.samstat_mapq_mqc.json
  ARGS --multiqc -d @OUTDIR@ ${SAMSTAT_FIXTURE_DIR}/small.sam
  REGEX "Failed to open output file")

# MultiQC file names are limited to 1024 characters: a -d path of 995 characters still
# takes the HTML report (+23) but not the longest MultiQC file (+35), which is an error
# and not a truncation. The directory is created below the scratch directory.
string(LENGTH "${MQC_OUT}/path_too_long" mqc_prefix_len)
math(EXPR mqc_rest "995 - ${mqc_prefix_len}")
set(MQC_LONG "")
while(mqc_rest GREATER 201)
  string(REPEAT "d" 200 comp)
  string(APPEND MQC_LONG "/${comp}")
  math(EXPR mqc_rest "${mqc_rest} - 201")
endwhile()
math(EXPR mqc_rest "${mqc_rest} - 1")
string(REPEAT "d" ${mqc_rest} comp)
string(APPEND MQC_LONG "/${comp}")
string(SUBSTRING "${MQC_LONG}" 1 -1 MQC_LONG_REL)
mqc_run_test(path_too_long 1
  MKDIR ${MQC_LONG_REL}
  ARGS --multiqc -d @OUTDIR@${MQC_LONG} ${SAMSTAT_FIXTURE_DIR}/small.sam
  REGEX "Output path too long")
