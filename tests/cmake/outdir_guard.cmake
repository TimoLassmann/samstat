# Included by every helper in this directory before it wipes OUTDIR.
# The helpers run file(REMOVE_RECURSE ${OUTDIR}) and let samstat write into OUTDIR, so a
# missing or wrong -DOUTDIR (a typo in a CMakeLists.txt, a variable that is not set) must
# stop the test instead of deleting or filling a directory outside the build tree.
# Every test directory is below SAMSTAT_TEST_OUT_DIR, <build>/tests/out (tests/CMakeLists.txt).

if(NOT DEFINED OUTDIR OR "${OUTDIR}" STREQUAL "")
  message(FATAL_ERROR "OUTDIR is not set: refusing to run (tests/README.md)")
endif()
if(NOT IS_ABSOLUTE "${OUTDIR}")
  message(FATAL_ERROR "OUTDIR '${OUTDIR}' is not an absolute path: refusing to run")
endif()
if(NOT "${OUTDIR}/" MATCHES "/tests/out/[^/]+/" OR "${OUTDIR}/" MATCHES "/\\.\\.?/")
  message(FATAL_ERROR "OUTDIR '${OUTDIR}' is not below <build>/tests/out: refusing to run")
endif()
