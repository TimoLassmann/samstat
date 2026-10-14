if(EXISTS "${FILE}")
  message(FATAL_ERROR "${FILE} must not exist in the source tree")
endif()
