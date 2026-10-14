# Run MultiQC on a directory of samstat _mqc files and check that every sample
# appears in each SAMStat section.
# Required: -DMULTIQC_EXE=<multiqc> -DINDIR=<dir> -DOUTDIR=<dir> -DSAMPLES="a;b;c"

file(REMOVE_RECURSE ${OUTDIR})
execute_process(
  COMMAND ${MULTIQC_EXE} ${INDIR} -o ${OUTDIR} -f
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "multiqc exited with ${rc}\nstdout:\n${out}\nstderr:\n${err}")
endif()

set(data ${OUTDIR}/multiqc_data)
set(files
  ${data}/multiqc_general_stats.txt
  ${data}/multiqc_samstat_mapq_table.txt
  ${data}/multiqc_samstat_length_plot.txt
  ${data}/multiqc_samstat_basequal_plot.txt)
foreach(f ${files})
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
endforeach()
