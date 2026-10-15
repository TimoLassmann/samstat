# Run MultiQC on a directory of samstat _mqc files and check that every sample
# appears in each SAMStat section.
# Required: -DMULTIQC_EXE=<multiqc> -DINDIR=<dir> -DOUTDIR=<dir> -DSAMPLES="a;b;c"
#
# MultiQC runs isolated from the user's setup: HOME, XDG_CONFIG_HOME and the
# matplotlib cache point into <OUTDIR>_env in the build tree, MULTIQC_CONFIG_PATH
# is unset, the working directory is <OUTDIR>_env (MultiQC reads a
# multiqc_config.yaml from the current directory) and a generated config turns
# off the version check, which would otherwise contact the MultiQC server.
# A multiqc installed with `pip install --user` keeps working on Linux and
# non-framework Pythons because PYTHONUSERBASE keeps the real ~/.local.

cmake_minimum_required(VERSION 3.19)

include(${CMAKE_CURRENT_LIST_DIR}/outdir_guard.cmake)
file(REMOVE_RECURSE ${OUTDIR})
set(envdir ${OUTDIR}_env)
file(REMOVE_RECURSE ${envdir})
file(MAKE_DIRECTORY ${envdir}/home ${envdir}/xdg ${envdir}/mpl)
set(config ${envdir}/multiqc_config.yaml)
file(WRITE ${config} "no_version_check: true\n")

# cmake -E env takes the --unset options before the NAME=VALUE pairs
set(env_args --unset=MULTIQC_CONFIG_PATH
  HOME=${envdir}/home XDG_CONFIG_HOME=${envdir}/xdg MPLCONFIGDIR=${envdir}/mpl)
if(NOT DEFINED ENV{PYTHONUSERBASE} AND DEFINED ENV{HOME})
  list(APPEND env_args PYTHONUSERBASE=$ENV{HOME}/.local)
endif()

execute_process(
  TIMEOUT 240 COMMAND ${CMAKE_COMMAND} -E env ${env_args}
    ${MULTIQC_EXE} -c ${config} ${INDIR} -o ${OUTDIR} -f
  WORKING_DIRECTORY ${envdir}
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
