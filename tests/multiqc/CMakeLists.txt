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

# Optional: a real MultiQC run. Registered only if multiqc is installed.
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
else()
  message(STATUS "multiqc not found: skipping optional multiqc:run test")
endif()
