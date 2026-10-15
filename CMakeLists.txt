cmake_minimum_required(VERSION 3.19)

project(SAMStat LANGUAGES C VERSION "2.4.0")

# C11 with the GNU/POSIX extensions (-std=gnu11): the sources use localtime_r, strdup/strndup,
# realpath, clock_gettime and getopt_long, which a strict -std=c11 does not declare (gcc 14 makes
# implicit declarations an error). Set here, not left to the default, so that a
# -DCMAKE_C_EXTENSIONS=OFF on the command line cannot turn them off.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # lots of warnings, but not (yet) as errors
  add_compile_options(-Wall -Wextra -pedantic)
  add_compile_options("$<$<CONFIG:DEBUG>:-O0>")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()


if(NOT PACKAGE_NAME)
set(PACKAGE_NAME ${CMAKE_PROJECT_NAME})
endif()

if(NOT PACKAGE_VERSION)
set(PACKAGE_VERSION ${CMAKE_PROJECT_VERSION})
endif()



option(SAMSTAT_SANITIZE "Build samstat with AddressSanitizer and UBSan" OFF)
option(SAMSTAT_USE_SYSTEM_HTSLIB "Use the system htslib (>=1.22) instead of the bundled one" OFF)

# config.h is generated into the build tree, never into the source tree.
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in"
            "${CMAKE_CURRENT_BINARY_DIR}/config.h" )
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# defines htslib::hts (bundled or system htslib)
include("thirdparty/htslib.cmake")

add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # Every test that sets no TIMEOUT of its own is killed after 120 s (the default is 1500 s):
  # a test that hangs or runs away must not keep running on a developer's machine.
  set(DART_TESTING_TIMEOUT 120 CACHE STRING "Default per-test timeout (seconds)" FORCE)
  include(CTest)
  # -DBUILD_TESTING=OFF (e.g. a package build) builds samstat only
  if(BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endif()

install(TARGETS samstat RUNTIME DESTINATION bin)

# Package Generator  #######################################################
# `cpack` makes a binary .tar.gz of the installed samstat; `cpack --config CPackSourceConfig.cmake`
# a source archive of the tree (with the bundled htslib tarball), without version control, build
# directories or editor files. Releases use the GitHub source archive (see release.yml).
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
# regular expressions on full paths; only top-level build directories are meant (tests/build is source)
string(REGEX REPLACE "([][+.*()^$?|\\\\])" "\\\\\\1" samstat_source_re "${CMAKE_CURRENT_SOURCE_DIR}")
string(REGEX REPLACE "([][+.*()^$?|\\\\])" "\\\\\\1" samstat_binary_re "${CMAKE_CURRENT_BINARY_DIR}")
# written to CPackSourceConfig.cmake as given, so each regex needs one level of escaping only
set(CPACK_VERBATIM_VARIABLES TRUE)
set(CPACK_SOURCE_IGNORE_FILES
  "/\\.git/"
  "/\\.git$"
  "/\\.cache/"
  "/\\.DS_Store$"
  "^${samstat_source_re}/build/"
  "^${samstat_source_re}/build-[^/]*/"
  "^${samstat_binary_re}/")
include(CPack)
