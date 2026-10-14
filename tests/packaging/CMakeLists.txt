add_test(NAME packaging:metadata
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPROJECT_VERSION=${PROJECT_VERSION}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/check_packaging.cmake)
set_tests_properties(packaging:metadata PROPERTIES LABELS packaging)
