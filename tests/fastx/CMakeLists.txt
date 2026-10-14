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
  "read1.length.unmapped.8=4")

# last record without trailing newline, gzip compressed (the plain variants exist above)
file(WRITE ${FASTX_GEN}/nonl_gz.fq "@r1\nACGT\n+\nIIII")
fastx_gzip(${FASTX_GEN}/nonl_gz.fq ${FASTX_GEN}/nonl.fq.gz)
fastx_path_test(fastx:json:nonl_fq_gz ${FASTX_GEN}/nonl.fq.gz
  "summary.total_reads=1"
  "read1.length.unmapped.4=1"
  "read1.base_quality.unmapped.40.3=1")
file(WRITE ${FASTX_GEN}/nonl_gz.fa ">a\nACGT\n>b\nACGT")
fastx_gzip(${FASTX_GEN}/nonl_gz.fa ${FASTX_GEN}/nonl.fa.gz)
fastx_path_test(fastx:json:nonl_fa_gz ${FASTX_GEN}/nonl.fa.gz
  "summary.total_reads=2"
  "read1.length.unmapped.4=2")

# FASTA sequences wrapped over several lines: s1 = 4 + 4 + 2 = 10 bp, s2 = 4 + 2 = 6 bp
file(WRITE ${FASTX_GEN}/multiline.fa ">s1 first sequence\nACGT\nACGT\nAC\n>s2\nGGGG\nCC\n")
fastx_path_test(fastx:json:multiline_fa ${FASTX_GEN}/multiline.fa
  "summary.total_reads=2"
  "summary.is_aligned=false"
  "read1.length.unmapped.10=1"
  "read1.length.unmapped.6=1")

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
  "read1.length.unmapped.12=1"
  "read1.length.unmapped.20=1")
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
  "read1.length.unmapped.4=2")

# white space only, compressed: an empty input
file(WRITE ${FASTX_GEN}/blank.txt "\n\n  \n")
fastx_gzip(${FASTX_GEN}/blank.txt ${FASTX_GEN}/blank.fq.gz)
fastx_path_test(fastx:blank_gz ${FASTX_GEN}/blank.fq.gz "summary.total_reads=0")
