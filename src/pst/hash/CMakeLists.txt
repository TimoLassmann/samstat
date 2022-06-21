
add_library(module_lhash STATIC
  lhash.c
  # lhash_io.c
  murmur3.c
  )

target_include_directories(module_lhash
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_lhash PRIVATE tld-dev)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# I could make this into a test... 
add_executable(unit_lhash
  unit_hash.c 
  )
target_link_libraries(unit_lhash PRIVATE tld-dev module_lhash)
