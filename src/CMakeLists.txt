add_executable(samstat
  samstat.c
  htsinterface/htsglue.c
  )

add_dependencies(samstat htslib)

target_link_libraries(samstat PRIVATE ${HTSLIB_SRC_DIR}/libhts.a tld-dev z bz2 curl lzma)

