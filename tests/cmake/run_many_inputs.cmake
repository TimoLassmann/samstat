# Several inputs at once: with at least as many inputs as threads each input runs on its own
# thread. The pages of -t 3 (4 inputs, 3 at a time) must be byte-identical to those of -t 1.
# The x-axis labels of plots 3-5 used to be static buffers shared by those threads, so an
# aligned input lost ", aligned part" whenever an unaligned one was drawn at the same moment.
# Two aligned (SAM) and two unaligned (FASTQ) inputs of 60 reads each, 10 rounds.
# Required: -DSAMSTAT=<binary> -DOUTDIR=<scratch dir, wiped first>

cmake_minimum_required(VERSION 3.19)

include(${CMAKE_CURRENT_LIST_DIR}/outdir_guard.cmake)
file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR}/in)

set(sam "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n")
set(fq "")
foreach(i RANGE 1 60)
  string(APPEND sam "r${i}\t0\tchr1\t${i}\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n")
  string(APPEND fq "@q${i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")
endforeach()
file(WRITE ${OUTDIR}/in/a.sam "${sam}")
file(WRITE ${OUTDIR}/in/b.sam "${sam}")
file(WRITE ${OUTDIR}/in/c.fq "${fq}")
file(WRITE ${OUTDIR}/in/d.fq "${fq}")
set(inputs ${OUTDIR}/in/a.sam ${OUTDIR}/in/b.sam ${OUTDIR}/in/c.fq ${OUTDIR}/in/d.fq)

file(MAKE_DIRECTORY ${OUTDIR}/ref)
execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} -t 1 -d ${OUTDIR}/ref ${inputs}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat -t 1 exited with ${rc}\n${out}\n${err}")
endif()
file(READ ${OUTDIR}/ref/a.sam.samstat.html ref_a)
string(FIND "${ref_a}" ", aligned part" at)
if(at EQUAL -1)
  message(FATAL_ERROR "-t 1: the position plots of a.sam do not say 'aligned part'; the test input is wrong")
endif()

foreach(round RANGE 1 10)
  file(MAKE_DIRECTORY ${OUTDIR}/r${round})
  execute_process(TIMEOUT 120 COMMAND ${SAMSTAT} -t 3 -d ${OUTDIR}/r${round} ${inputs}
    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "samstat -t 3 exited with ${rc}\n${out}\n${err}")
  endif()
  foreach(name a.sam b.sam c.fq d.fq)
    file(READ ${OUTDIR}/ref/${name}.samstat.html want)
    file(READ ${OUTDIR}/r${round}/${name}.samstat.html got)
    if(NOT want STREQUAL got)
      message(FATAL_ERROR "round ${round}: ${name}.samstat.html differs between -t 1 and -t 3")
    endif()
  endforeach()
  file(REMOVE_RECURSE ${OUTDIR}/r${round})
endforeach()
file(REMOVE_RECURSE ${OUTDIR}/in)
