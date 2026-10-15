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
# With offset 0 the '!' would be Q33 and the 'I' Q73. E[q] = round(10^(-q/10) * 2^30):
# E[0] = 1073741824, E[40] = 107374. Position 0 sums E[0] + E[40], positions 1..3 2 * E[40].
# Per-read mean Q: r2 is Q40; r1 has mean error (E[0] + 3 E[40]) / 4 = 268515986, the largest
# q with E[q] >= 268515986 is 6.
fastx_json_test(qual qual.fq
  "summary.total_reads=2"
  "summary.is_aligned=false"
  "read1.absolute.qual_bases.unmapped#len=4"
  "read1.absolute.qual_bases.unmapped.0=2"
  "read1.absolute.qual_bases.unmapped.3=2"
  "read1.absolute.qual_error_sum.unmapped.0=1073849198"
  "read1.absolute.qual_error_sum.unmapped.1=214748"
  "read1.absolute.qual_error_sum.unmapped.3=214748"
  "read1.read_mean_q.unmapped#len=41"
  "read1.read_mean_q.unmapped.6=1"
  "read1.read_mean_q.unmapped.40=1")

# inputs that are valid but contain no sequence
fastx_json_test(empty_gz empty.fq.gz "summary.total_reads=0")
fastx_json_test(newline_only newline_only.fq "summary.total_reads=0")
fastx_json_test(header_only header_only.fa "summary.total_reads=0")

# the last record has no trailing newline
fastx_json_test(nonl_fq nonl.fq
  "summary.total_reads=1"
  "read1.length.reads.unmapped.4=1"
  "read1.absolute.qual_bases.unmapped.3=1"
  "read1.absolute.qual_error_sum.unmapped.3=107374")
fastx_json_test(nonl_fa nonl.fa
  "summary.total_reads=2"
  "read1.length.reads.unmapped.4=2")

add_test(NAME fastx:warn:empty_gz
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DARGS=-d\;@OUTDIR@\;${SAMSTAT_FIXTURE_DIR}/empty.fq.gz"
    -DOUTDIR=${FASTX_OUT}/warn_empty_gz "-DEXPECT_STDERR_REGEX=no reads"
    "-DEXPECT_FILES=empty.fq.gz.samstat.html"
    -P ${FASTX_HELPERS}/run_expect.cmake)
set_tests_properties(fastx:warn:empty_gz PROPERTIES LABELS fastx)

# =============================================================================
# WP2.1: edge cases and pending tests. Inputs that are not in tests/data are written
# into the build tree at configure time (nothing is committed for them).
# =============================================================================
set(FASTX_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${FASTX_GEN})

# fastx_path_test(<test name> <input path> <checks...>): value checks on the JSON report
function(fastx_path_test name input)
  string(REPLACE ";" "\;" checks "${ARGN}")
  string(REPLACE ":" "_" dir "${name}")
  add_test(NAME ${name}
    COMMAND ${CMAKE_COMMAND}
      -DSAMSTAT=$<TARGET_FILE:samstat> -DINPUT=${input}
      -DOUTDIR=${FASTX_OUT}/${dir} "-DCHECKS=${checks}"
      -P ${FASTX_HELPERS}/run_json.cmake)
  set_tests_properties(${name} PROPERTIES LABELS fastx)
endfunction()

# gzip a file at configure time (CMake's own libarchive, no external gzip needed)
function(fastx_gzip src dst)
  file(ARCHIVE_CREATE OUTPUT ${dst} PATHS ${src} FORMAT raw COMPRESSION GZip)
endfunction()

# --- enabled now (pass on 2.3) -----------------------------------------------------------
# Read names ending /1 and /2 do not make FASTQ reads paired: every read is read 1. Pins
# the 2.3 behaviour; the htslib reader (WP3.2) sets BAM_FREAD2 from the suffix and must
# ignore it for FASTQ.
file(WRITE ${FASTX_GEN}/r1r2.fq
  "@pairA/1\nACGTACGT\n+\nIIIIIIII\n@pairA/2\nTTGGCCAA\n+\nIIIIIIII\n"
  "@pairB/1\nACGTACGT\n+\nIIIIIIII\n@pairB/2\nTTGGCCAA\n+\nIIIIIIII\n")
fastx_path_test(fastx:r1r2_suffix ${FASTX_GEN}/r1r2.fq
  "summary.total_reads=4"
  "summary.read1_reads=4"
  "summary.read2_reads=0"
  "read1.length.reads.unmapped.8=4")

# last record without trailing newline, gzip compressed (the plain variants exist above)
file(WRITE ${FASTX_GEN}/nonl_gz.fq "@r1\nACGT\n+\nIIII")
fastx_gzip(${FASTX_GEN}/nonl_gz.fq ${FASTX_GEN}/nonl.fq.gz)
fastx_path_test(fastx:json:nonl_fq_gz ${FASTX_GEN}/nonl.fq.gz
  "summary.total_reads=1"
  "read1.length.reads.unmapped.4=1"
  "read1.absolute.qual_bases.unmapped.3=1"
  "read1.absolute.qual_error_sum.unmapped.3=107374")
file(WRITE ${FASTX_GEN}/nonl_gz.fa ">a\nACGT\n>b\nACGT")
fastx_gzip(${FASTX_GEN}/nonl_gz.fa ${FASTX_GEN}/nonl.fa.gz)
fastx_path_test(fastx:json:nonl_fa_gz ${FASTX_GEN}/nonl.fa.gz
  "summary.total_reads=2"
  "read1.length.reads.unmapped.4=2")

# FASTA sequences wrapped over several lines: s1 = 4 + 4 + 2 = 10 bp, s2 = 4 + 2 = 6 bp
file(WRITE ${FASTX_GEN}/multiline.fa ">s1 first sequence\nACGT\nACGT\nAC\n>s2\nGGGG\nCC\n")
fastx_path_test(fastx:json:multiline_fa ${FASTX_GEN}/multiline.fa
  "summary.total_reads=2"
  "summary.is_aligned=false"
  "read1.length.reads.unmapped.10=1"
  "read1.length.reads.unmapped.6=1")

# --- content detection (the type comes from the content, WP3.2) --------------------------
# The type comes from the content, not from the extension. Copies are made at configure
# time (macOS file systems are case-insensitive: the names differ by more than case).
configure_file(${SAMSTAT_FIXTURE_DIR}/READS.FQ ${FASTX_GEN}/misnamed_fastq.sam COPYONLY)
configure_file(${SAMSTAT_FIXTURE_DIR}/small.sam ${FASTX_GEN}/misnamed_sam.fq COPYONLY)
configure_file(${SAMSTAT_FIXTURE_DIR}/reads.fq.gz ${FASTX_GEN}/misnamed_gz.txt COPYONLY)
fastx_path_test(fastx:sniff_misnamed ${FASTX_GEN}/misnamed_fastq.sam
  "summary.total_reads=3"
  "summary.is_aligned=false")
fastx_path_test(fastx:sniff_misnamed_sam ${FASTX_GEN}/misnamed_sam.fq
  "summary.total_reads=20"
  "summary.is_aligned=true")
fastx_path_test(fastx:sniff_misnamed_gz ${FASTX_GEN}/misnamed_gz.txt
  "summary.total_reads=3"
  "summary.is_aligned=false")

# Protein FASTA is unsupported. htslib takes a file for FASTA only if its first sequence
# line is made of IUPAC nucleotide letters; such a file is read as nucleotides (letters
# other than ACGT count as N). Any other protein FASTA is rejected: exit 1 and a message.
file(WRITE ${FASTX_GEN}/iupac_protein.fa ">p1\nMKRSACGTNNRY\n>p2\nMKRSACGTNNRYMKRSACGT\n")
fastx_path_test(fastx:protein_fasta_iupac ${FASTX_GEN}/iupac_protein.fa
  "summary.total_reads=2"
  "summary.is_aligned=false"
  "read1.length.reads.unmapped.12=1"
  "read1.length.reads.unmapped.20=1")
file(WRITE ${FASTX_GEN}/protein.fa ">p1\nMKVLAAGIVGLLLAQ\n>p2\nMSTNPKPQRKTKRNTNRRPQ\n")
add_test(NAME fastx:protein_fasta_rejected
  COMMAND ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DARGS=-d\;@OUTDIR@\;${FASTX_GEN}/protein.fa"
    -DOUTDIR=${FASTX_OUT}/protein_fasta_rejected -DEXPECT_RC=1 "-DEXPECT_STDERR_REGEX=protein FASTA is not supported"
    "-DEXPECT_NO_FILES=protein.fa.samstat.html"
    -P ${FASTX_HELPERS}/run_expect.cmake)
set_tests_properties(fastx:protein_fasta_rejected PROPERTIES LABELS fastx)

# a record with an empty sequence is not a read
file(WRITE ${FASTX_GEN}/empty_seq.fa ">r0\nACGT\n>e1\n>r1\nACGT\n>e2\n")
fastx_path_test(fastx:empty_sequence_fa ${FASTX_GEN}/empty_seq.fa
  "summary.total_reads=2"
  "read1.length.reads.unmapped.4=2")

# white space only, compressed: an empty input
file(WRITE ${FASTX_GEN}/blank.txt "\n\n  \n")
fastx_gzip(${FASTX_GEN}/blank.txt ${FASTX_GEN}/blank.fq.gz)
fastx_path_test(fastx:blank_gz ${FASTX_GEN}/blank.fq.gz "summary.total_reads=0")

# --- streaming pre-check (precheck.c): a record or line too long is rejected before htslib
# loads it. SAMSTAT_MAX_READ_LEN lowers the read length limit to 30 so that the inputs stay
# tiny; a line may be 65536 characters longer than that (PRECHECK_LINE_SLACK). The message
# differs from the reader's own check ("has N bases, more than ..."), so these tests fail if
# the input reaches htslib's record reader first.
# fastx_limit_test(<name> <input> RC <n> [REGEX <re>] [CHECKS <json checks...>])
function(fastx_limit_test name input)
  cmake_parse_arguments(T "" "RC;REGEX" "CHECKS" ${ARGN})
  get_filename_component(base ${input} NAME)
  string(REPLACE ";" "\;" checks "${T_CHECKS}")
  set(cmd ${CMAKE_COMMAND} -E env SAMSTAT_MAX_READ_LEN=30 ${CMAKE_COMMAND}
    -DSAMSTAT=$<TARGET_FILE:samstat> "-DARGS=--json\;-d\;@OUTDIR@\;${input}"
    -DOUTDIR=${FASTX_OUT}/precheck_${name} -DEXPECT_RC=${T_RC})
  if(T_REGEX)
    list(APPEND cmd "-DEXPECT_STDERR_REGEX=${T_REGEX}")
  endif()
  if(T_RC EQUAL 0)
    list(APPEND cmd "-DEXPECT_FILES=${base}.samstat.html" -DJSON=${base}.samstat.json "-DCHECKS=${checks}")
  else()
    list(APPEND cmd "-DEXPECT_NO_FILES=${base}.samstat.html")
  endif()
  add_test(NAME fastx:precheck:${name} COMMAND ${cmd} -P ${FASTX_HELPERS}/run_expect.cmake)
  set_tests_properties(fastx:precheck:${name} PROPERTIES LABELS fastx)
endfunction()

# FASTA: the second record has 20 + 20 = 40 bases, over the limit of 30 (plain and gzip)
file(WRITE ${FASTX_GEN}/pc_big.fa
  ">ok\nACGTACGTAC\n>big chromosome\nACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\n")
fastx_gzip(${FASTX_GEN}/pc_big.fa ${FASTX_GEN}/pc_big.fa.gz)
fastx_limit_test(fasta_record ${FASTX_GEN}/pc_big.fa RC 1
  REGEX "record big is longer than 30 bases: not a file of sequencing reads")
fastx_limit_test(fasta_record_gz ${FASTX_GEN}/pc_big.fa.gz RC 1
  REGEX "record big is longer than 30 bases: not a file of sequencing reads")

# FASTQ: one record of 2 x 20 = 40 bases on two sequence and two quality lines
file(WRITE ${FASTX_GEN}/pc_big.fq
  "@q1\nACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\nIIIIIIIIIIIIIIIIIIII\n")
fastx_limit_test(fastq_multiline ${FASTX_GEN}/pc_big.fq RC 1
  REGEX "record q1 is longer than 30 bases: not a file of sequencing reads")

# gzip: one line of 70001 characters and no newline (htslib takes it for FASTQ and would
# read the line whole); over 30 + 65536
string(REPEAT "x" 70000 pc_blob)
file(WRITE ${FASTX_GEN}/pc_long_line.txt "@${pc_blob}")
fastx_gzip(${FASTX_GEN}/pc_long_line.txt ${FASTX_GEN}/pc_long_line.fq.gz)
fastx_limit_test(gz_long_line ${FASTX_GEN}/pc_long_line.fq.gz RC 1
  REGEX "line 1 is longer than 65566 characters: not a file of sequencing reads")

# SAM: a header line over the line limit
file(WRITE ${FASTX_GEN}/pc_long_header.sam "@HD\tVN:1.6\n@CO\t${pc_blob}\n")
fastx_limit_test(sam_header_line ${FASTX_GEN}/pc_long_header.sam RC 1
  REGEX "line 2 is longer than 65566 characters: not a file of sequencing reads")

# accepted: records of exactly 30 bases. FASTQ: sequence and quality over two lines, the
# quality lines start with '@' and '+' (they are quality, not a new record: htslib reads
# quality lines until they cover the sequence), a name line with a 200-character comment
# (longer than 30, within the slack). FASTA: CRLF line ends (the '\r' is not a base).
string(REPEAT "c" 200 pc_comment)
file(WRITE ${FASTX_GEN}/pc_ok.fq
  "@r1 ${pc_comment}\nACGTACGTACGTACG\nTACGTACGTACGTAC\n+\n@IIIIIIIIIIIIII\n+IIIIIIIIIIIIII\n"
  "@r2\nACGTACGTACGTACGTACGTACGTACGTAC\n+\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
  "@r3\nACGT\n+\nIIII\n")
fastx_gzip(${FASTX_GEN}/pc_ok.fq ${FASTX_GEN}/pc_ok.fq.gz)
fastx_limit_test(accept_fq ${FASTX_GEN}/pc_ok.fq RC 0 CHECKS
  "summary.total_reads=3"
  "read1.length.reads.unmapped.30=2"
  "read1.length.reads.unmapped.4=1")
fastx_limit_test(accept_fq_gz ${FASTX_GEN}/pc_ok.fq.gz RC 0 CHECKS
  "summary.total_reads=3"
  "read1.length.reads.unmapped.30=2"
  "read1.length.reads.unmapped.4=1")
file(WRITE ${FASTX_GEN}/pc_ok.fa
  ">s1 ${pc_comment}\r\nACGTACGTACGTACG\r\nTACGTACGTACGTAC\r\n>s2\r\nACGTACGTACGTACGTACGTACGTACGTAC\r\n")
fastx_limit_test(accept_fa_crlf ${FASTX_GEN}/pc_ok.fa RC 0 CHECKS
  "summary.total_reads=2"
  "read1.length.reads.unmapped.30=2")
