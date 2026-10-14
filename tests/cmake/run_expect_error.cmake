# Run samstat on an input that must be rejected and check the diagnostics.
# Required: -DSAMSTAT=<binary> -DINPUT=<file> -DOUTDIR=<dir>
#           -DEXPECT_REGEX=<regex that must appear in stdout/stderr>
# The test fails if samstat crashes (signal, "Segmentation", "Abort"...) or if a
# report is written for the rejected input. The exit code is not checked here:
# it is covered by the exit-code tests.

cmake_minimum_required(VERSION 3.19)

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

execute_process(
  TIMEOUT 120 COMMAND ${SAMSTAT} --json -d ${OUTDIR} ${INPUT}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)

# a crash is reported by execute_process as a non-numeric result such as
# "Segmentation fault" or "Child aborted"
if(NOT rc MATCHES "^[0-9]+$")
  message(FATAL_ERROR "samstat crashed on ${INPUT}: ${rc}\nstdout:\n${out}\nstderr:\n${err}")
endif()

if(NOT "${out}${err}" MATCHES "${EXPECT_REGEX}")
  message(FATAL_ERROR "Expected message matching '${EXPECT_REGEX}' for ${INPUT}\nexit: ${rc}\nstdout:\n${out}\nstderr:\n${err}")
endif()

file(GLOB reports ${OUTDIR}/*.html ${OUTDIR}/*.json)
if(reports)
  message(FATAL_ERROR "Report written for rejected input ${INPUT}: ${reports}")
endif()
