add_test(NAME packaging:metadata
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPROJECT_VERSION=${PROJECT_VERSION}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/check_packaging.cmake)
set_tests_properties(packaging:metadata PROPERTIES LABELS packaging)

add_test(NAME packaging:tld_pin
  COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/check_tld_pin.cmake)
set_tests_properties(packaging:tld_pin PROPERTIES LABELS packaging)
