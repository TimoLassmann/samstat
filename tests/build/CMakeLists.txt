# Build-system tests (label: build)

add_test(NAME build:version
  COMMAND samstat -v)
set_tests_properties(build:version PROPERTIES
  LABELS build
  PASS_REGULAR_EXPRESSION "${PROJECT_VERSION}")

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
