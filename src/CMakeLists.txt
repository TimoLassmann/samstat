add_executable(samstat
  samstat.c
  htsinterface/htsglue.c
  sambamparse/sam_bam_parse.c
  )


target_link_directories(samstat PUBLIC ${HTSLIB_SRC_DIR})


target_link_libraries(samstat PRIVATE z)
target_link_libraries(samstat PRIVATE bz2)
target_link_libraries(samstat PRIVATE lzma)
target_link_libraries(samstat PRIVATE curl)

if(APPLE)
  target_include_directories(samstat PRIVATE ${HTSLIB_SRC_DIR})
  target_link_libraries(samstat PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)
else()
  target_link_libraries(samstat PRIVATE hts)
endif()


#target_link_libraries(samstat PRIVATE hts)
target_link_libraries(samstat PRIVATE tld-dev)
