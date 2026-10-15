# Build-system tests (label: build)

string(REPLACE "." "\\." PROJECT_VERSION_RE "${PROJECT_VERSION}")
add_test(NAME build:version
  COMMAND samstat -v)
set_tests_properties(build:version PROPERTIES
  LABELS build
  PASS_REGULAR_EXPRESSION "${PROJECT_VERSION_RE}")

add_test(NAME build:help
  COMMAND samstat -h)
set_tests_properties(build:help PROPERTIES
  LABELS build
  PASS_REGULAR_EXPRESSION "Usage")

add_test(NAME build:no_config_h_in_source
  COMMAND ${CMAKE_COMMAND} -DFILE=${PROJECT_SOURCE_DIR}/src/config.h
    -P ${CMAKE_CURRENT_SOURCE_DIR}/no_file.cmake)
set_tests_properties(build:no_config_h_in_source PROPERTIES LABELS build)

add_test(NAME build:config_h_in_build_tree
  COMMAND ${CMAKE_COMMAND} -DFILE=${PROJECT_BINARY_DIR}/config.h
    -P ${CMAKE_CURRENT_SOURCE_DIR}/has_file.cmake)
set_tests_properties(build:config_h_in_build_tree PROPERTIES LABELS build)

add_test(NAME build:no_hdf5_link
  COMMAND ${CMAKE_COMMAND} -DSAMSTAT=$<TARGET_FILE:samstat>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/no_hdf5_link.cmake)
set_tests_properties(build:no_hdf5_link PROPERTIES LABELS build)

# build:convert_tables_in_sync  the committed src/collect/convert_tables.h is what
# src/convert_tables.c generates. The generator is built for the tests only (samstat no longer
# runs it) and has to run on the build host, so the test is skipped in a cross build.
if(NOT CMAKE_CROSSCOMPILING)
  add_executable(converttablegen ${PROJECT_SOURCE_DIR}/src/convert_tables.c)
  target_link_libraries(converttablegen PRIVATE tldlite)
  add_test(NAME build:convert_tables_in_sync
    COMMAND ${CMAKE_COMMAND}
      -DGEN=$<TARGET_FILE:converttablegen>
      -DOUT=${CMAKE_CURRENT_BINARY_DIR}/convert_tables.h
      -DREF=${PROJECT_SOURCE_DIR}/src/collect/convert_tables.h
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tables_in_sync.cmake)
  set_tests_properties(build:convert_tables_in_sync PROPERTIES LABELS build TIMEOUT 120)
endif()
