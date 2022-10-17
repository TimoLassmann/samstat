
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endif()

# if(NOT EXISTS "${PROJECT_SOURCE_DIR}/thirdparty/htslib/Makefile")
#   message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
# endif()

add_subdirectory(tld)


add_executable(samstat
  samstat.c
  htsinterface/htsglue.c
  sambamparse/sam_bam_parse.c
  param/param.c
  metrics/metrics.c
  report/report.c
  report/lst.c
  tools/tools.c
  )

target_link_directories(samstat PUBLIC ${HTSLIB_SRC_DIR})


# target_link_libraries(samstat PRIVATE z)
# target_link_libraries(samstat PRIVATE bz2)
# target_link_libraries(samstat PRIVATE lzma)
# target_link_libraries(samstat PRIVATE curl)

if(APPLE)
  target_include_directories(samstat PRIVATE ${HTSLIB_SRC_DIR})
  target_link_libraries(samstat PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)
else()
  target_include_directories(samstat PRIVATE ${HTSLIB_SRC_DIR})
  target_link_libraries(samstat PRIVATE hts)
endif()


#target_link_libraries(samstat PRIVATE hts)
target_link_libraries(samstat PRIVATE tld-dev)
