include(ExternalProject)


set(config_flags)  # parameters desired for ./configure of Autotools

#set(my_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}my${CMAKE_STATIC_LIBRARY_SUFFIX})

find_program(MAKE_EXECUTABLE NAMES gmake make mingw32-make REQUIRED)

# set(HTSLIB_SOURCES_DIR ${CMAKE_CURRENT_LIST_DIR}/htslib)
# set(HTSLIB_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/htslib/htslib)

set(htslib_PREFIX ${CMAKE_SOURCE_DIR}/thirdparty/htslib)
include_directories(${htslib_PREFIX}/htslib)
set(htslib_LIB ${htslib_PREFIX}/libhts.a)

ExternalProject_Add(htslib
  GIT_REPOSITORY https://github.com/samtools/htslib.git
  GIT_TAG 1.15.1
  GIT_SUBMODULES htscodecs

  BUILD_IN_SOURCE 1
  BUILD_ALWAYS 1
  #UPDATE_DISCONNECTED true
  #CONFIGURE_HANDLED_BY_BUILD false
  #BUILD_IN_SOURCE 1
  SOURCE_DIR "${htslib_PREFIX}"
  CONFIGURE_COMMAND ${htslib_PREFIX}/configure ${config_flags}
  BUILD_COMMAND ${MAKE_EXECUTABLE} -j
  INSTALL_COMMAND ""  
  TEST_COMMAND ""
  #BUILD_BYPRODUCTS 
  )

ExternalProject_Add_Step(
  htslib autoheader
  COMMAND "autoheader"
  WORKING_DIRECTORY "${htslib_PREFIX}"
  DEPENDERS configure
  DEPENDEES autoconf
  )

ExternalProject_Add_Step(
  htslib autoconf
  COMMAND "autoconf"  
  WORKING_DIRECTORY "${htslib_PREFIX}"
  DEPENDERS autoheader
  DEPENDEES download
  )

# ExternalProject_Get_Property(htslib source_dir)
# set(HTSLIB_INCLUDE_DIR ${source_dir}/htslib)
# include_directories(${HTSLIB_INCLUDE_DIR})
# ExternalProject_Get_Property(htslib binary_dir)
# add_library(libhts STATIC IMPORTED)
# set(LIBHTS_A ${binary_dir}/libhts.a)
# set_property(TARGET libhts PROPERTY IMPORTED_LOCATION ${LIBHTS_A})

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

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/thirdparty/htslib/Makefile")
  message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()


add_subdirectory(tld)
