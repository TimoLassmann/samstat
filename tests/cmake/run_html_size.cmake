# The HTML report stays below 1 MB in both layouts (standard-plots.md section 1).
# Writes the input into OUTDIR, runs samstat and checks the size of the report.
# Required: -DSAMSTAT=<binary> -DOUTDIR=<dir>
#           -DMODE=short       a paired, aligned SAM of 999 bp reads with enough records to
#                              give every position bin of every MAPQ class and mate data
#                              (worst case for the number of plotted points, short layout)
#           -DMODE=ultralong   one aligned read of 1 Mb (long layout)
# Optional: -DMAX_BYTES=<n>    default 1000000

cmake_minimum_required(VERSION 3.19)

if(NOT DEFINED MAX_BYTES)
  set(MAX_BYTES 1000000)
endif()

file(REMOVE_RECURSE ${OUTDIR})
file(MAKE_DIRECTORY ${OUTDIR})

set(header "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:2000000\n")

if(MODE STREQUAL "ultralong")
  set(len 1000000)
  math(EXPR quarter "${len}/4")
  string(REPEAT "ACGT" ${quarter} seq)
  string(REPEAT "I" ${len} qual)
  set(sam ${OUTDIR}/ultralong.sam)
  file(WRITE ${sam} "${header}"
    "ultra\t0\tchr1\t1\t60\t${len}M\t*\t0\t0\t${seq}\t${qual}\tNM:i:0\tMD:Z:${len}\n")
elseif(MODE STREQUAL "short")
  # Pair i: MAPQ class i%4 (60, 20, 5, unmapped), one mismatch (MD) at read position i,
  # base qualities that differ per pair, so that every position bin of every class has
  # data. 1000 pairs cover all 999 positions.
  set(len 999)
  string(REPEAT "C" ${len} seq)
  set(sam ${OUTDIR}/short_999.sam)
  file(WRITE ${sam} "${header}")
  set(mapqs 60 20 5 0)
  set(chunk "")
  foreach(i RANGE 0 999)
    math(EXPR cls "${i}%4")
    math(EXPR pos "1+(${i}*37)%90000")
    math(EXPR mate "${pos}+200")
    math(EXPR mm "${i}%${len}")
    math(EXPR rest "${len}-1-${mm}")
    math(EXPR q "33+(${i}*7)%41")
    string(ASCII ${q} qc)
    math(EXPR q2 "33+(${i}*11+5)%41")
    string(ASCII ${q2} qc2)
    string(REPEAT "${qc}" ${len} qual1)
    string(REPEAT "${qc2}" ${len} qual2)
    if(cls EQUAL 3)
      # unmapped pair
      string(APPEND chunk
        "p${i}\t77\t*\t0\t0\t*\t*\t0\t0\t${seq}\t${qual1}\n"
        "p${i}\t141\t*\t0\t0\t*\t*\t0\t0\t${seq}\t${qual2}\n")
    else()
      list(GET mapqs ${cls} mapq)
      string(APPEND chunk
        "p${i}\t99\tchr1\t${pos}\t${mapq}\t${len}M\t=\t${mate}\t1200\t${seq}\t${qual1}\tNM:i:1\tMD:Z:${mm}A${rest}\n"
        "p${i}\t147\tchr1\t${mate}\t${mapq}\t${len}M\t=\t${pos}\t-1200\t${seq}\t${qual2}\tNM:i:1\tMD:Z:${mm}A${rest}\n")
    endif()
    math(EXPR m "${i}%100")
    if(m EQUAL 99)
      file(APPEND ${sam} "${chunk}")
      set(chunk "")
    endif()
  endforeach()
  file(APPEND ${sam} "${chunk}")
else()
  message(FATAL_ERROR "run_html_size: set -DMODE=short or -DMODE=ultralong")
endif()

execute_process(COMMAND ${SAMSTAT} -d ${OUTDIR} ${sam}
  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "samstat exited with ${rc} (${MODE})\nstdout:\n${out}\nstderr:\n${err}")
endif()
get_filename_component(name ${sam} NAME)
set(report ${OUTDIR}/${name}.samstat.html)
if(NOT EXISTS ${report})
  message(FATAL_ERROR "No HTML report written: ${report}")
endif()
file(SIZE ${report} size)
if(NOT size LESS ${MAX_BYTES})
  message(FATAL_ERROR "${MODE}: HTML report is ${size} bytes, expected < ${MAX_BYTES}")
endif()
message(STATUS "${MODE}: HTML report ${size} bytes (< ${MAX_BYTES})")
file(REMOVE ${sam})
