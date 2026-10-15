# samstat must not link HDF5 (the tld library once pulled it in; tldlite never does).
# Required: -DSAMSTAT=<binary>
if(APPLE)
  set(cmd otool -L ${SAMSTAT})
else()
  set(cmd ldd ${SAMSTAT})
endif()
execute_process(COMMAND ${cmd} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err
  TIMEOUT 120)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${cmd} failed (${rc}):\n${err}")
endif()
string(TOLOWER "${out}" out_lc)
if(out_lc MATCHES "hdf5")
  message(FATAL_ERROR "samstat links HDF5:\n${out}")
endif()
