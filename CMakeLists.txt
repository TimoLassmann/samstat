cmake_minimum_required(VERSION 3.2)

project(SAMStat LANGUAGES C VERSION "2.2.1")

# message(STATUS "${CMAKE_PROJECT_NAME} Name ")
# message(STATUS "${CMAKE_PROJECT_VERSION} Version ")


set (CMAKE_C_STANDARD 11)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

add_compile_options("$<$<CONFIG:RELEASE>:-W;-Wall;-pedantic>")
add_compile_options("$<$<CONFIG:DEBUG>:-W;-Wall;-O0;-pedantic;-DDEBUG>")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if (MSVC)
    # warning level 4 and all warnings as errors
    add_compile_options(/W4 /WX)
  else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic) # -Werror )
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#include(GNUInstallDirs)

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



configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/config.h" )
include(ExternalProject)
include("thirdparty/htslib.cmake")

#option(BUILD_SHARED_LIBS "Build all libraries as shared" OFF)
add_subdirectory(thirdparty)


add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)  
#  add_subdirectory(tests)
endif()
#add_subdirectory(src/tld/tests)

install(TARGETS samstat RUNTIME DESTINATION bin)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
include(CPack)  
# Package Generator  #######################################################
# set(CPACK_PACKAGE_VENDOR "Vendor name")
# set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Some summary")
# set(CPACK_PACKAGE_CONTACT "Timo Lassmann")
# set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
# set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
# set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
# set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
# set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
# set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
# set(CPACK_SOURCE_IGNORE_FILES
#      /.git
#      /dist
#      /GPATH
#      /GTAGS
#      /GRTAGS
#      /.*build.*
#      /.dir-locals.el
#      /\\\\.DS_Store
#  )
#  include (CPack)


