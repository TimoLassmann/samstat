
add_executable(samstat
  samstat.c
  )





set(htslib_PREFIX ${CMAKE_SOURCE_DIR}/thirdparty/htslib/htslib)

set(htslib_LIB ${CMAKE_SOURCE_DIR}/thirdparty/htslib/libhts.so)

# if(NOT EXISTS ${htslib_PREFIX}/sam.h)
#   message(SEND_ERROR "Can't find sam.h in ${htslib_PREFIX}")
# endif()

#include_directories(${htslib_PREFIX})



#target_include_directories(samstat PUBLIC ${htslib_PREFIX}/htslib)


# add_dependencies(samstat htslib)
target_include_directories(samstat PRIVATE ${htslib_PREFIX}) 
message(STATUS "INCLUDE::::  ${htslib_PREFIX}")

link_directories(${CMAKE_SOURCE_DIR}/thirdparty/htslib)

find_library(FOO_LIB hts )

message(STATUS "INCLUDE::::  ${FOO_LIB}")
target_link_libraries(samstat PRIVATE tld-dev "${FOO_LIB}" z m bz2 lzma pthread curl)



#target_link_libraries(samstat PRIVATE ${HTSLIB} )  # need the quotes to expand

