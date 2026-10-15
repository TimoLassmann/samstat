# Command line, threading and exit code tests (label: cli).
#
# Exit codes: 0 success (also for empty inputs, with a warning), 1 at least one
# input failed (the other inputs are still processed), 2 usage error.
# All tests use tests/cmake/run_expect.cmake.

set(CLI_HELPERS ${PROJECT_SOURCE_DIR}/tests/cmake)

# the test-only environment limits can only lower the built-in values (param.c)
add_executable(unit_param unit_param.c ${PROJECT_SOURCE_DIR}/src/param/param.c)
target_include_directories(unit_param PRIVATE ${PROJECT_BINARY_DIR})
target_link_libraries(unit_param PRIVATE tldlite)
add_test(NAME cli:unit_param COMMAND unit_param)
set_tests_properties(cli:unit_param PROPERTIES LABELS cli)
set(CLI_OUT ${SAMSTAT_TEST_OUT_DIR}/cli)
set(CLI_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${CLI_GEN})

# cli_test(<name> ARGS <args...> [RC n] [REGEX re] [NO_REGEX re] [FILES f...]
#          [NO_FILES f...] [GEN_SAM n] [JSON file CHECKS c...] [NO_NAN])
# Arguments may use @OUTDIR@ (scratch dir of the test) and @GEN@ (generated SAM) and @GENFQ@ (generated FASTQ).
function(cli_test name)
  cmake_parse_arguments(T "NO_NAN" "RC;REGEX;NO_REGEX;GEN_SAM;GEN_FASTQ;GEN_BLOB;JSON" "ARGS;FILES;NO_FILES;CHECKS;ENV" ${ARGN})
  if(NOT DEFINED T_RC)
    set(T_RC 0)
  endif()
  string(REPLACE ";" "\;" args "${T_ARGS}")
  string(REPLACE ";" "\;" files "${T_FILES}")
  string(REPLACE ";" "\;" nofiles "${T_NO_FILES}")
  string(REPLACE ";" "\;" checks "${T_CHECKS}")
  # ENV: test-only variables for samstat (they can only lower built-in limits, see param.c)
  set(cmd ${CMAKE_COMMAND})
  if(T_ENV)
    list(APPEND cmd -E env ${T_ENV} ${CMAKE_COMMAND})
  endif()
  list(APPEND cmd
    -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${CLI_OUT}/${name}
    "-DARGS=${args}" -DEXPECT_RC=${T_RC})
  if(T_REGEX)
    list(APPEND cmd "-DEXPECT_STDERR_REGEX=${T_REGEX}")
  endif()
  if(T_NO_REGEX)
    list(APPEND cmd "-DEXPECT_NO_REGEX=${T_NO_REGEX}")
  endif()
  if(T_FILES)
    list(APPEND cmd "-DEXPECT_FILES=${files}")
  endif()
  if(T_NO_FILES)
    list(APPEND cmd "-DEXPECT_NO_FILES=${nofiles}")
  endif()
  if(T_GEN_SAM)
    list(APPEND cmd -DGEN_SAM=${T_GEN_SAM})
  endif()
  if(T_GEN_FASTQ)
    list(APPEND cmd -DGEN_FASTQ=${T_GEN_FASTQ})
  endif()
  if(T_GEN_BLOB)
    list(APPEND cmd -DGEN_BLOB=${T_GEN_BLOB})
  endif()
  if(T_JSON)
    list(APPEND cmd "-DJSON=${T_JSON}" "-DCHECKS=${checks}")
  endif()
  if(T_NO_NAN)
    list(APPEND cmd -DNO_NAN=ON)
  endif()
  list(APPEND cmd -P ${CLI_HELPERS}/run_expect.cmake)
  add_test(NAME cli:${name} COMMAND ${cmd})
  set_tests_properties(cli:${name} PROPERTIES LABELS cli)
endfunction()

set(FX ${SAMSTAT_FIXTURE_DIR})
set(D @OUTDIR@)

# --- worker error handling / exit codes of failing inputs ---------------------
# 1 = input failed, other inputs still processed
foreach(threads 1 4)
  cli_test(mixed_good_bad_t${threads}
    ARGS -t ${threads} -d ${D} ${FX}/small.sam ${FX}/does_not_exist.sam
    RC 1 REGEX "does_not_exist\\.sam"
    FILES small.sam.samstat.html)
  cli_test(mixed_bad_good_t${threads}
    ARGS -t ${threads} -d ${D} ${FX}/does_not_exist.sam ${FX}/small.sam
    RC 1 REGEX "does_not_exist\\.sam"
    FILES small.sam.samstat.html)
endforeach()
cli_test(nonexistent_file
  ARGS -d ${D} ${FX}/does_not_exist.bam
  RC 1 REGEX "does_not_exist\\.bam")
cli_test(all_good_rc0
  ARGS -t 2 -d ${D} ${FX}/small.sam ${FX}/small.bam
  FILES small.sam.samstat.html small.bam.samstat.html)

# --- usage errors exit 2, help / version exit 0 --------------------------------
cli_test(no_args ARGS "" RC 2 REGEX "Usage")
cli_test(help_short ARGS -h RC 0 REGEX "Usage")
cli_test(help_long ARGS -help RC 0 REGEX "Usage")
cli_test(version_short ARGS -v RC 0 REGEX "[0-9]+\\.[0-9]+")
cli_test(version_long ARGS -version RC 0 REGEX "[0-9]+\\.[0-9]+")
cli_test(threads_long ARGS --threads 2 -d ${D} ${FX}/small.sam FILES small.sam.samstat.html)
cli_test(nthreads_alias ARGS --nthreads 2 -d ${D} ${FX}/small.sam FILES small.sam.samstat.html)
cli_test(only_options ARGS -t 2 RC 2 REGEX "Usage")
cli_test(unknown_option ARGS --bogus ${FX}/small.sam RC 2 NO_FILES small.sam.samstat.html)
cli_test(missing_option_value ARGS ${FX}/small.sam -t RC 2)
cli_test(duplicate_input
  ARGS -d ${D} ${FX}/small.sam ${FX}/small.sam
  RC 2 REGEX "Duplicated input" NO_FILES small.sam.samstat.html)
# -d creates the directory (mkdir -p)
cli_test(dir_created
  ARGS -d ${D}/new/deeper ${FX}/small.sam
  FILES new/deeper/small.sam.samstat.html)

# --- file type detection --------------------------------------------------------
# The type comes from the content, never from the extension (WP3.2). An input whose
# content nothing recognises fails only that input: exit 1, a message naming it, and
# the other inputs still get their report.
file(WRITE ${CLI_GEN}/garbage.xyz "this is not a sequence file\nat all\n")
cli_test(ext_unknown_only ARGS -d ${D} ${CLI_GEN}/garbage.xyz RC 1 REGEX "garbage\\.xyz"
  NO_FILES garbage.xyz.samstat.html)
cli_test(ext_unknown_with_good
  ARGS -d ${D} ${FX}/small.sam ${CLI_GEN}/garbage.xyz
  RC 1 REGEX "garbage\\.xyz" FILES small.sam.samstat.html NO_FILES garbage.xyz.samstat.html)
cli_test(ext_fq_gz ARGS -d ${D} ${FX}/reads.fq.gz FILES reads.fq.gz.samstat.html)
cli_test(ext_fa_gz ARGS -d ${D} ${FX}/reads.fa.gz FILES reads.fa.gz.samstat.html)
cli_test(ext_upper_case ARGS -d ${D} ${FX}/READS.FQ FILES READS.FQ.samstat.html)
# a file that does not exist fails that input (exit 1) whatever its name
cli_test(ext_not_a_suffix_only
  ARGS -d ${D} ${FX}/small.sam.txt RC 1 REGEX "small\\.sam\\.txt")

# upper/mixed case and the .gz variants of every fasta/fastq extension, via copies
set(EXT_COPY_DIR ${CLI_OUT}/ext_copies)
file(MAKE_DIRECTORY ${EXT_COPY_DIR})
# NOTE: macOS file systems are case-insensitive, so every test needs an output
# directory name that differs by more than case (hence the counter).
set(n_ext 0)
foreach(ext SAM Sam BAM Bam)
  math(EXPR n_ext "${n_ext}+1")
  if(ext MATCHES "^[Ss][Aa][Mm]$")
    set(src ${FX}/small.sam)
  else()
    set(src ${FX}/small.bam)
  endif()
  configure_file(${src} ${EXT_COPY_DIR}/case.${ext} COPYONLY)
  cli_test(ext_case_${n_ext}_${ext} ARGS -d ${D} ${EXT_COPY_DIR}/case.${ext} FILES case.${ext}.samstat.html)
endforeach()
foreach(ext fa fasta fna fq fastq FA FASTQ)
  math(EXPR n_ext "${n_ext}+1")
  configure_file(${FX}/READS.FQ ${EXT_COPY_DIR}/plain.${ext} COPYONLY)
  cli_test(ext_plain_${n_ext}_${ext} ARGS -d ${D} ${EXT_COPY_DIR}/plain.${ext} FILES plain.${ext}.samstat.html)
endforeach()
foreach(ext fa.gz fasta.gz fna.gz fq.gz fastq.gz FQ.GZ Fa.Gz)
  math(EXPR n_ext "${n_ext}+1")
  configure_file(${FX}/reads.fq.gz ${EXT_COPY_DIR}/gz.${ext} COPYONLY)
  cli_test(ext_gz_${n_ext}_${ext} ARGS -d ${D} ${EXT_COPY_DIR}/gz.${ext} FILES gz.${ext}.samstat.html)
endforeach()

# --- numeric options are validated strictly (exit 2) --------------------------------
foreach(opt -t -p)
  foreach(bad abc 12x -5 0x10 " 7" 1.5 99999999999999999999999)
    string(MAKE_C_IDENTIFIER "${opt}_${bad}" id)
    cli_test(badnum_${id}
      ARGS ${opt} "${bad}" ${FX}/small.sam
      RC 2 REGEX "Option" NO_FILES small.sam.samstat.html)
  endforeach()
endforeach()
cli_test(badnum_t_zero ARGS -t 0 ${FX}/small.sam RC 2 REGEX "-t")
cli_test(badnum_p_zero ARGS -p 0 ${FX}/small.sam RC 2 REGEX "-p")
# -t is capped at 64: larger values are rejected before any thread or buffer is created.
# Never test with a large accepted-looking value (e.g. -t 100000): if the cap ever broke, the
# test itself would start that many threads. 65 is just over the cap; a value beyond 32 bits
# is rejected by the number parser, independently of the cap.
cli_test(threads_max_ok ARGS -d ${D} -t 64 ${FX}/small.sam FILES small.sam.samstat.html)
cli_test(threads_over_max_65 ARGS -d ${D} -t 65 ${FX}/small.sam
  RC 2 REGEX "at most 64 threads" NO_FILES small.sam.samstat.html)
cli_test(threads_not_a_32bit_number ARGS -d ${D} -t 99999999999 ${FX}/small.sam
  RC 2 REGEX "too large" NO_FILES small.sam.samstat.html)

# --- resource limits (memory must not grow with a setting or a malformed input) ---------------
# a record longer than the read length limit fails its input (built-in 50 Mb; lowered here so that
# the test input stays tiny)
# (limit 30: the 20 bp reads of the generated FASTQ pass, the 40 bp reads of small.sam fail)
cli_test(read_longer_than_limit ENV SAMSTAT_MAX_READ_LEN=30
  GEN_FASTQ 3 ARGS -d ${D} @GENFQ@ ${FX}/small.sam
  RC 1 REGEX "record small00 has 40 bases, more than the 30 a read may have"
  FILES gen_3.fq.samstat.html NO_FILES small.sam.samstat.html)
cli_test(read_at_limit ENV SAMSTAT_MAX_READ_LEN=20
  GEN_FASTQ 3 ARGS -d ${D} @GENFQ@ FILES gen_3.fq.samstat.html)
# (that the test-only limits can only lower the built-in values is checked by cli:unit_param,
# which calls param_init() directly and never allocates by the values)
# 1 MB without a newline that is not SAM/FASTA/FASTQ: rejected after reading a fixed chunk,
# never loaded as one line
cli_test(blob_without_newline GEN_BLOB 1000000 ARGS -d ${D} @GENBLOB@
  RC 1 REGEX "unrecognised input")
cli_test(num_ok ARGS -d ${D} -t 3 -p 100 ${FX}/small.sam FILES small.sam.samstat.html)

# --- -p/--peek analyses exactly N reads (reader batch is 4096 records) --------------------------
# peek_test(<name> <kind> <n_reads> <peek> <expected total> <partial true|false>)
function(peek_test name kind n_reads peek expected partial)
  if(kind STREQUAL "sam")
    set(gen GEN_SAM ${n_reads})
    set(file @GEN@)
    set(out gen_${n_reads}.sam.samstat.json)
  else()
    set(gen GEN_FASTQ ${n_reads})
    set(file @GENFQ@)
    set(out gen_${n_reads}.fq.samstat.json)
  endif()
  if(partial)
    set(rx "first ${expected} reads")
  else()
    set(rx "")
  endif()
  cli_test(peek_${name} ${gen} ARGS --json -d ${D} -p ${peek} ${file}
    JSON ${out} CHECKS summary.total_reads=${expected} summary.is_partial_report=${partial}
    REGEX "${rx}")
endfunction()
peek_test(sam_10_of_30   sam 30 10 10 true)
peek_test(sam_1_of_30    sam 30 1 1 true)
peek_test(sam_all        sam 30 30 30 false)
peek_test(sam_more_than  sam 30 1000 30 false)
peek_test(sam_batch_edge sam 4100 4096 4096 true)
peek_test(sam_batch_plus sam 4100 4097 4097 true)
peek_test(sam_batch_all  sam 4100 4100 4100 false)
peek_test(sam_chunk_edge sam 10050 10000 10000 true)
peek_test(sam_chunk_plus sam 10050 10010 10010 true)
peek_test(sam_exact_two_chunks sam 10050 10050 10050 false)
peek_test(fq_10_of_30    fq 30 10 10 true)
peek_test(fq_all         fq 30 30 30 false)
peek_test(fq_batch_edge  fq 4100 4096 4096 true)
peek_test(fq_batch_plus  fq 4100 4097 4097 true)
peek_test(fq_batch_all   fq 4100 4100 4100 false)
peek_test(fq_chunk_edge  fq 10050 10000 10000 true)
peek_test(fq_chunk_plus  fq 10050 10010 10010 true)
peek_test(fq_exact_two_chunks fq 10050 10050 10050 false)

# --- empty inputs: report written, warning, exit 0, no nan/inf in the output ----------
foreach(f empty.sam)
  cli_test(empty_sam
    ARGS --json -d ${D} ${FX}/${f}
    REGEX "[Ee]mpty|no reads|0 reads" NO_NAN
    FILES empty.sam.samstat.html empty.sam.samstat.json
    JSON empty.sam.samstat.json CHECKS summary.total_reads=0)
endforeach()
cli_test(empty_fq
  ARGS --json -d ${D} ${FX}/empty.fq
  REGEX "[Ee]mpty|no reads|0 reads" NO_NAN
  FILES empty.fq.samstat.html empty.fq.samstat.json
  JSON empty.fq.samstat.json CHECKS summary.total_reads=0)
cli_test(empty_with_good_t1
  ARGS -t 1 -d ${D} ${FX}/empty.sam ${FX}/small.sam
  FILES empty.sam.samstat.html small.sam.samstat.html)

# --- error paths free their resources and report failure ----------------------------
# (run under AddressSanitizer / LeakSanitizer in CI to catch leaks on these paths)
cli_test(error_fastq_no_sequence
  ARGS -d ${D} ${CMAKE_CURRENT_SOURCE_DIR}/data/bad_no_seq.fq
  RC 1 REGEX "bad_no_seq\\.fq")
cli_test(error_fastq_truncated
  ARGS -d ${D} ${CMAKE_CURRENT_SOURCE_DIR}/data/bad_truncated.fq ${FX}/small.sam
  RC 1 REGEX "bad_truncated\\.fq")
cli_test(error_good_after_bad_fastq
  ARGS -t 1 -d ${D} ${CMAKE_CURRENT_SOURCE_DIR}/data/bad_no_seq.fq ${FX}/small.sam
  RC 1 FILES small.sam.samstat.html)

# --- threads give the same results ---------------------------------------------
add_test(NAME cli:threads_many_inputs
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${FX}/small.sam
    -DOUTDIR=${CLI_OUT}/threads_many_inputs -P ${CLI_HELPERS}/run_threads.cmake)
set_tests_properties(cli:threads_many_inputs PROPERTIES LABELS cli)

# --- WP1.4: deprecated options, up-front output checks, version and help -------------
# The scripts are generated into the build tree so that tests/cmake/ stays untouched.

# deprecated option: exit 0, exactly one warning line, JSON identical to a run without it
file(WRITE ${CLI_GEN}/run_deprecated.cmake [=[
# -DSAMSTAT -DOUTDIR -DINPUT -DOPT=<name> -DOPTARG=<value or empty>
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR}/a ${OUTDIR}/b)
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${OUTDIR}/a ${INPUT}
  RESULT_VARIABLE rc_a OUTPUT_VARIABLE out_a ERROR_VARIABLE err_a)
set(flag ${OPT})
if(OPTARG)
  list(APPEND flag ${OPTARG})
endif()
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} ${flag} --json -d ${OUTDIR}/b ${INPUT}
  RESULT_VARIABLE rc_b OUTPUT_VARIABLE out_b ERROR_VARIABLE err_b)
if(NOT rc_a EQUAL 0 OR NOT rc_b EQUAL 0)
  message(FATAL_ERROR "exit codes ${rc_a} / ${rc_b}\n${out_b}\n${err_b}")
endif()
set(log "${out_b}\n${err_b}")
string(REGEX MATCHALL "samstat: warning: [^\n]*deprecated[^\n]*" w "${log}")
list(LENGTH w nw)
if(NOT nw EQUAL 1)
  message(FATAL_ERROR "expected exactly one deprecation warning, got ${nw}:\n${log}")
endif()
if(NOT w MATCHES "samstat: warning: ${OPT} is deprecated and ignored")
  message(FATAL_ERROR "unexpected warning text: ${w}")
endif()
file(READ ${OUTDIR}/a/small.sam.samstat.json ja)
file(READ ${OUTDIR}/b/small.sam.samstat.json jb)
if(NOT ja STREQUAL jb)
  message(FATAL_ERROR "JSON differs when ${OPT} is given")
endif()
]=])
function(deprecated_test name opt optarg)
  add_test(NAME cli:deprecated_${name}
    COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${CLI_OUT}/deprecated_${name}
      -DINPUT=${FX}/small.sam -DOPT=${opt} "-DOPTARG=${optarg}" -P ${CLI_GEN}/run_deprecated.cmake)
  set_tests_properties(cli:deprecated_${name} PROPERTIES LABELS cli)
endfunction()
deprecated_test(l -l 100)
deprecated_test(seed --seed 7)
deprecated_test(norm --norm "")
deprecated_test(plotend --plotend "")

# unwritable report directories fail that input (exit 1), quickly, naming -d DIR
file(WRITE ${CLI_GEN}/run_perm.cmake [=[
# -DSAMSTAT -DOUTDIR -DINPUT -DMODE=outdir|input_dir
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR}/ro)
if(MODE STREQUAL "outdir")
  set(args -d ${OUTDIR}/ro ${INPUT})
else()
  file(COPY_FILE ${INPUT} ${OUTDIR}/ro/small.sam)
  set(args ${OUTDIR}/ro/small.sam)
endif()
file(CHMOD ${OUTDIR}/ro PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
# root ignores permission bits: skip
execute_process(COMMAND id -u OUTPUT_VARIABLE uid OUTPUT_STRIP_TRAILING_WHITESPACE)
if(uid STREQUAL "0")
  file(CHMOD ${OUTDIR}/ro PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  file(REMOVE_RECURSE ${OUTDIR})
  message(STATUS "running as root, skipping")
  return()
endif()
string(TIMESTAMP t0 "%s")
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} ${args} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
string(TIMESTAMP t1 "%s")
file(CHMOD ${OUTDIR}/ro PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
math(EXPR dt "${t1}-${t0}")
set(msg "")
if(NOT rc EQUAL 1)
  set(msg "exit code ${rc}, expected 1")
elseif(NOT "${out}${err}" MATCHES "cannot write [^\n]*; use -d DIR")
  set(msg "message lacks 'use -d DIR'")
elseif(dt GREATER 1)
  set(msg "took ${dt} s")
endif()
if(MODE STREQUAL "outdir")
  file(GLOB made ${OUTDIR}/ro/*.samstat.html)
  if(made)
    set(msg "${msg} report written")
  endif()
endif()
# leave nothing behind that a backup tool or ctest cleanup would trip over
file(REMOVE_RECURSE ${OUTDIR})
if(msg)
  message(FATAL_ERROR "${msg}\nstdout:\n${out}\nstderr:\n${err}")
endif()
]=])
foreach(mode outdir input_dir)
  set(tname outdir_unwritable)
  if(mode STREQUAL "input_dir")
    set(tname input_dir_readonly)
  endif()
  add_test(NAME cli:${tname}
    COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${CLI_OUT}/perm_${mode}
      -DINPUT=${FX}/small.sam -DMODE=${mode} -P ${CLI_GEN}/run_perm.cmake)
  set_tests_properties(cli:${tname} PROPERTIES LABELS cli)
endforeach()

# two inputs resolving to the same report file: exit 2, both named
file(MAKE_DIRECTORY ${CLI_OUT}/collide/a ${CLI_OUT}/collide/b)
file(COPY_FILE ${FX}/small.sam ${CLI_OUT}/collide/a/same.sam ONLY_IF_DIFFERENT)
file(COPY_FILE ${FX}/small.sam ${CLI_OUT}/collide/b/same.sam ONLY_IF_DIFFERENT)
cli_test(output_collision
  ARGS -d ${D} ${CLI_OUT}/collide/a/same.sam ${CLI_OUT}/collide/b/same.sam
  RC 2 REGEX "collide/a/same\\.sam.*collide/b/same\\.sam" NO_FILES same.sam.samstat.html)

cli_test(version_regex ARGS -v RC 0 REGEX "(^|\n)SAMStat [0-9]+\\.[0-9]+\\.[0-9]+\n")

# every long option of -h is documented in Readme.org
add_test(NAME cli:help_matches_readme
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DREADME=${PROJECT_SOURCE_DIR}/Readme.org
    -P ${CLI_GEN}/check_help.cmake)
set_tests_properties(cli:help_matches_readme PROPERTIES LABELS cli)
file(WRITE ${CLI_GEN}/check_help.cmake [=[
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} -h OUTPUT_VARIABLE help RESULT_VARIABLE rc)
file(READ ${README} readme)
string(REGEX MATCHALL "--[a-z]+" opts "${help}")
if(NOT opts)
  message(FATAL_ERROR "no long options found in -h output")
endif()
foreach(o IN LISTS opts)
  string(FIND "${readme}" "${o}" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "${o} is in -h but not in Readme.org")
  endif()
endforeach()
# the usage line and the option lines are the same text in both, in the same order
string(REGEX MATCH "Usage: [a-z]+ ([^\n]*)" usage "${help}")
string(FIND "${readme}" "samstat ${CMAKE_MATCH_1}\n" at)
if(at EQUAL -1)
  message(FATAL_ERROR "the usage line '${usage}' of -h is not in Readme.org")
endif()
string(REGEX MATCHALL "\n +-[^\n]*" lines "${help}")
set(from 0)
foreach(l IN LISTS lines)
  string(SUBSTRING "${l}" 1 -1 l)
  string(FIND "${readme}" "${l}\n" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "the -h line '${l}' is not in Readme.org")
  endif()
  if(at LESS from)
    message(FATAL_ERROR "the -h line '${l}' is in a different order in Readme.org")
  endif()
  set(from ${at})
endforeach()
]=])

# =============================================================================
# Threading, -p and standard input tests.
# =============================================================================

# --- WP4.1 (intra-file threading): one file, -t 1/2/3/4/8 give identical output ------------
# JSON, HTML and the four MultiQC files must be byte-identical; --verbose reports the thread
# budget ('workers: W, decompression threads: P'). SAMSTAT_BATCH_RECORDS=8 and
# SAMSTAT_BATCH_BASES=20000 (test-only) cut the input into at least 5 batches.
# paired_small.sam has only 4 records, so its records are repeated 40 times (160 records,
# 20 batches). ont_like.bam: 25 reads, 89 kb, 5 or more batches (BAM, decompression pool).
# long.fq: FASTQ (no alignment data). reads.fq.gz: gzip FASTQ through the pool.
add_test(NAME cli:threads_one_file:paired_small
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${FX}/paired_small.sam
    -DOUTDIR=${CLI_OUT}/threads_one_file_paired -DMODE=one_file -DREPLICATE=40
    -P ${CLI_HELPERS}/run_threads.cmake)
add_test(NAME cli:threads_one_file:ont_like
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${FX}/longread/ont_like.bam
    -DOUTDIR=${CLI_OUT}/threads_one_file_ont -DMODE=one_file
    -P ${CLI_HELPERS}/run_threads.cmake)
add_test(NAME cli:threads_one_file:long_fq
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${FX}/longread/long.fq
    -DOUTDIR=${CLI_OUT}/threads_one_file_long_fq -DMODE=one_file
    -P ${CLI_HELPERS}/run_threads.cmake)
add_test(NAME cli:threads_one_file:reads_fq_gz
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${FX}/reads.fq.gz
    -DOUTDIR=${CLI_OUT}/threads_one_file_fq_gz -DMODE=one_file
    -P ${CLI_HELPERS}/run_threads.cmake)
set_tests_properties(cli:threads_one_file:paired_small cli:threads_one_file:ont_like
  cli:threads_one_file:long_fq cli:threads_one_file:reads_fq_gz
  PROPERTIES LABELS cli)

# The -p cutoff at -t 8: the same cases and input as reader:p_cutoff_counts (tests/reader,
# which owns the generated script), with the reader and 5 workers on one file.
set(P_CUTOFF_CASES_T8
  "1:1:0:0:true" "6:6:1:1:true" "8:8:1:1:true" "4096:4096:512:512:true"
  "4097:4097:512:512:true" "5003:5003:626:625:true" "7199:7199:900:900:true"
  "7200:7200:900:900:false" "100000:7200:900:900:false")
string(REPLACE ";" "\;" p_cutoff_cases_t8 "${P_CUTOFF_CASES_T8}")
add_test(NAME reader:p_cutoff_counts_t8
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${CLI_OUT}/p_cutoff_counts_t8
    -DCHECK_JSON=${CLI_HELPERS}/check_json.cmake -DTHREADS=8 "-DCASES=${p_cutoff_cases_t8}"
    -P ${PROJECT_BINARY_DIR}/tests/reader/gen/run_p_cutoff.cmake)
set_tests_properties(reader:p_cutoff_counts_t8 PROPERTIES LABELS reader)

# --- WP3.2 (htslib reader, content detection, stdin) -----------------------------------
# stdin: the input is fed through standard input with `-` as the file name. The reports are
# named stdin.samstat.*; most tests glob for the one JSON report in -d DIR and assert on its
# content, cli:stdin_report_names pins the names.
file(WRITE ${CLI_GEN}/run_stdin.cmake [=[
# -DSAMSTAT -DOUTDIR -DSTDIN_FILE=<file fed to stdin> -DCHECKS=<check_json list>
# -DCHECK_JSON=<path of check_json.cmake>
# optional: -DEXPECT_RC=<n> (default 0); with a nonzero code only the exit code is checked
#           -DEXTRA_ARGS=<more arguments after the first "-">
#           -DEXPECT_FILES=<list of files in OUTDIR that must exist>
#           -DEXPECT_ERR_REGEX=<regex the log output (stdout and stderr) must match>
#           -DHTML_TEXT=<text stdin.samstat.html must contain>
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})
if(NOT DEFINED EXPECT_RC)
  set(EXPECT_RC 0)
endif()
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${OUTDIR} - ${EXTRA_ARGS}
  INPUT_FILE ${STDIN_FILE}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT "${rc}" STREQUAL "${EXPECT_RC}")
  message(FATAL_ERROR "samstat - < ${STDIN_FILE}: exit ${rc}, expected ${EXPECT_RC}\nstdout:\n${out}\nstderr:\n${err}")
endif()
if(NOT EXPECT_RC EQUAL 0)
  return()
endif()
string(REPLACE "\\;" ";" EXPECT_FILES "${EXPECT_FILES}")
foreach(f IN LISTS EXPECT_FILES)
  if(NOT EXISTS ${OUTDIR}/${f})
    message(FATAL_ERROR "expected file ${f} was not written to ${OUTDIR}\nstdout:\n${out}\nstderr:\n${err}")
  endif()
endforeach()
if(EXPECT_ERR_REGEX AND NOT "${out}\n${err}" MATCHES "${EXPECT_ERR_REGEX}")
  message(FATAL_ERROR "log output does not match '${EXPECT_ERR_REGEX}'\nstdout:\n${out}\nstderr:\n${err}")
endif()
if(HTML_TEXT)
  file(READ ${OUTDIR}/stdin.samstat.html html)
  string(FIND "${html}" "${HTML_TEXT}" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "stdin.samstat.html does not contain '${HTML_TEXT}'")
  endif()
endif()
file(GLOB reports ${OUTDIR}/*.json)
list(LENGTH reports n)
if(NOT n EQUAL 1)
  message(FATAL_ERROR "expected exactly one JSON report in ${OUTDIR}, found ${n}: ${reports}\nstdout:\n${out}\nstderr:\n${err}")
endif()
execute_process(
  COMMAND ${CMAKE_COMMAND} -DJSON=${reports} "-DCHECKS=${CHECKS}"
    -P ${CHECK_JSON}
  RESULT_VARIABLE rc2 OUTPUT_VARIABLE out2 ERROR_VARIABLE err2)
if(NOT rc2 EQUAL 0)
  message(FATAL_ERROR "${out2}${err2}")
endif()
]=])
# stdin_test(<name> <fixture> <checks...>)
function(stdin_test name input)
  string(REPLACE ";" "\;" checks "${ARGN}")
  add_test(NAME cli:stdin_${name}
    COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
      -DOUTDIR=${CLI_OUT}/stdin_${name} -DSTDIN_FILE=${input} -DCHECK_JSON=${CLI_HELPERS}/check_json.cmake "-DCHECKS=${checks}"
      -P ${CLI_GEN}/run_stdin.cmake)
  set_tests_properties(cli:stdin_${name} PROPERTIES LABELS cli)
endfunction()
stdin_test(sam ${FX}/small.sam "summary.total_reads=20" "summary.is_aligned=true")
stdin_test(bam ${FX}/small.bam "summary.total_reads=20" "summary.is_aligned=true")
stdin_test(fastq ${FX}/READS.FQ "summary.total_reads=3" "summary.is_aligned=false")
stdin_test(fastq_gz ${FX}/reads.fq.gz "summary.is_aligned=false" "summary.total_reads=3")
stdin_test(paired_sam ${FX}/paired_small.sam
  "summary.total_reads=4" "summary.read1_reads=2" "summary.read2_reads=2")

# the report names for standard input (decision C1), the stderr summary line, and the format
# in the page header: stdin is a stream, the format comes from the reader's own detection
add_test(NAME cli:stdin_report_names
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${CLI_OUT}/stdin_report_names -DSTDIN_FILE=${FX}/small.bam -DCHECK_JSON=${CLI_HELPERS}/check_json.cmake
    "-DCHECKS=summary.total_reads=20" "-DEXPECT_FILES=stdin.samstat.html\;stdin.samstat.json"
    "-DEXPECT_ERR_REGEX=samstat: stdin -> [^\n]*stdin\\.samstat\\.html \\(BAM, 20 reads"
    "-DHTML_TEXT=>BAM &middot; " -P ${CLI_GEN}/run_stdin.cmake)
set_tests_properties(cli:stdin_report_names PROPERTIES LABELS cli)

# empty standard input is a valid empty input (3.0): exit 0, total_reads 0
add_test(NAME cli:stdin_empty
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${CLI_OUT}/stdin_empty -DSTDIN_FILE=${FX}/zero_bytes.sam -DCHECK_JSON=${CLI_HELPERS}/check_json.cmake
    "-DCHECKS=summary.total_reads=0" -P ${CLI_GEN}/run_stdin.cmake)
set_tests_properties(cli:stdin_empty PROPERTIES LABELS cli)

# stdin with content nothing recognises fails that input: exit 1, no report
file(WRITE ${CLI_GEN}/garbage.txt "this is not a sequence file\nat all\n")
add_test(NAME cli:stdin_garbage
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${CLI_OUT}/stdin_garbage -DSTDIN_FILE=${CLI_GEN}/garbage.txt
    -DEXPECT_RC=1 -P ${CLI_GEN}/run_stdin.cmake)
set_tests_properties(cli:stdin_garbage PROPERTIES LABELS cli)

# an input whose content is not recognised fails only that input (exit 1), whatever its
# extension; the other inputs still get their report. (unknown.xyz is a valid SAM, so it
# is accepted by content: see cli:unknown_input_sniffed.)
cli_test(unknown_input_per_file
  ARGS -d ${D} ${FX}/small.sam ${CLI_GEN}/garbage.xyz
  RC 1 REGEX "garbage\\.xyz"
  FILES small.sam.samstat.html NO_FILES garbage.xyz.samstat.html)
cli_test(unknown_input_sniffed
  ARGS -d ${D} ${FX}/unknown.xyz
  FILES unknown.xyz.samstat.html)

# standard input may be given only once: usage error, exit 2, nothing is read
add_test(NAME cli:stdin_double_dash
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -DOUTDIR=${CLI_OUT}/stdin_double_dash -DSTDIN_FILE=${FX}/small.sam
    -DEXPECT_RC=2 -DEXTRA_ARGS=- -P ${CLI_GEN}/run_stdin.cmake)
set_tests_properties(cli:stdin_double_dash PROPERTIES LABELS cli)

# --- stderr summary line and progress (standard-plots.md section 6) ---------------------
# one line per input: samstat: <in> -> <out> (<format>, <n> reads, N50 ..., mean Q ..., [% mapped,] <layout> layout, <t> s)
cli_test(stderr_summary_line_bam
  ARGS -d ${D} ${FX}/small.bam
  REGEX "samstat: [^\n]*small\\.bam -> [^\n]*small\\.bam\\.samstat\\.html \\(BAM, [0-9.]+[kMG]? reads, N50 [0-9.]+ (bp|kb|Mb), mean Q [0-9.]+, [0-9.]+% mapped, (short|long) layout, [0-9]+ s\\)")
cli_test(stderr_summary_line_fastq
  GEN_FASTQ 50 ARGS -d ${D} @GENFQ@
  REGEX "samstat: [^\n]*gen_50\\.fq -> [^\n]*gen_50\\.fq\\.samstat\\.html \\(FASTQ, 50 reads, N50 20 bp, mean Q 40\\.0, short layout, [0-9]+ s\\)"
  NO_REGEX "% mapped")
cli_test(stderr_summary_line_partial
  GEN_FASTQ 50 ARGS -d ${D} -p 10 @GENFQ@
  REGEX "FASTQ, 10 reads, N50 20 bp, mean Q 40\\.0, short layout, partial: first 10 reads, [0-9]+ s\\)")
# a pipe is not a terminal: no counter ("N reads processed" + CR)
# (the progress step is lowered to 100 reads, so a small file crosses it many times)
cli_test(no_progress_without_tty ENV SAMSTAT_PROGRESS_STEP=100
  GEN_FASTQ 1100 ARGS -d ${D} @GENFQ@
  REGEX "FASTQ, 1\\.1k reads"
  NO_REGEX "reads processed|Processed")

# --- final review fixes ------------------------------------------------------------------

# Output files are not announced on stdout unless --verbose is given (standard-plots.md section 6)
cli_test(no_writing_to_log
  ARGS --json --multiqc -d ${D} ${FX}/small.sam
  FILES small.sam.samstat.html small.sam.samstat.json small.sam.samstat_general_mqc.json
  NO_REGEX "Writing to")
cli_test(verbose_writing_to_log
  ARGS --verbose --json --multiqc -d ${D} ${FX}/small.sam
  REGEX "Writing to: [^\n]*small\\.sam\\.samstat\\.html.*Writing to: [^\n]*small\\.sam\\.samstat\\.json.*Writing to: [^\n]*samstat_general_mqc\\.json")

# JSON mapq_bins labels read as the key of the MAPQ table does
cli_test(json_mapq_bin_labels
  ARGS --json -d ${D} ${FX}/small.sam
  JSON small.sam.samstat.json
  CHECKS "mapq_bins.0.label=MAPQ >= 30" "mapq_bins.1.label=MAPQ 10-29" "mapq_bins.2.label=MAPQ < 10" "mapq_bins.3.label=Unmapped")

# A mapped record with SEQ '*' has no bases: its soft clips must not enter the clipped-bases
# percentage (it was 1000 % for 100 clipped bases over 10 mapped ones)
file(WRITE ${CLI_GEN}/clip_seq_star.sam
  "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n"
  "normal\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
  "star\t0\tchr1\t20\t60\t100S10M\t*\t0\t0\t*\t*\n")
cli_test(pct_clipped_seq_star
  ARGS --json -d ${D} ${CLI_GEN}/clip_seq_star.sam
  JSON clip_seq_star.sam.samstat.json CHECKS "summary.pct_clipped<=100" "summary.pct_clipped>=0")

# -p: empty-sequence FASTQ records are not reads, so they do not make the report partial
file(WRITE ${CLI_GEN}/empty_seq_last.fq "@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n@c\n\n+\n\n")
cli_test(peek_empty_seq_not_partial
  ARGS --json -p 2 -d ${D} ${CLI_GEN}/empty_seq_last.fq
  NO_REGEX "partial: |partial report"
  JSON empty_seq_last.fq.samstat.json CHECKS "summary.total_reads=2" "summary.is_partial_report=false")
cli_test(peek_real_record_is_partial
  ARGS --json -p 1 -d ${D} ${CLI_GEN}/empty_seq_last.fq
  REGEX "partial: first 1 reads"
  JSON empty_seq_last.fq.samstat.json CHECKS "summary.total_reads=1" "summary.is_partial_report=true")

# Several inputs processed at once (-t 3, 4 inputs): the HTML is byte-identical to -t 1. The axis
# labels of plots 3-5 were once shared between threads.
add_test(NAME cli:threads_many_inputs_html
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat> -DOUTDIR=${CLI_OUT}/threads_many_inputs_html
    -P ${CLI_HELPERS}/run_many_inputs.cmake)
set_tests_properties(cli:threads_many_inputs_html PROPERTIES LABELS cli)
