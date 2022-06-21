
add_subdirectory(hash)


 add_library(module_pst STATIC
   pst.c
   pst_alloc.c
  )

target_include_directories(module_pst
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_pst PRIVATE tld-dev module_lhash)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# I could make this into a test... 
# add_executable(unit_pst
#   pst_model_test.c
#   )
# target_link_libraries(unit_pst PRIVATE tld-dev module_pst_core)
