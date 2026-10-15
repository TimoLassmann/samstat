add_subdirectory(util)
add_subdirectory(plot)

option(USE_PTHREADS    "Use pthreads for parallelization" ON)

if(USE_PTHREADS)
  find_package(Threads REQUIRED)
  if(Threads_FOUND)
    message(STATUS "pthreads is enabled.")
    add_definitions (-DHAVE_PTHREADS)
  else(Threads_FOUND)
    message(STATUS "pthreads not enabled")
  endif(Threads_FOUND)
endif(USE_PTHREADS)

add_executable(converttablegen
  convert_tables.c
  )
target_link_libraries(converttablegen PRIVATE tldlite)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/convert_tables.h
  COMMAND converttablegen ${CMAKE_CURRENT_BINARY_DIR}/convert_tables.h
  DEPENDS converttablegen
  )

set(source_files
  samstat.c
  htsinterface/htsglue.c
  htsinterface/precheck.c
  sambamparse/sam_bam_parse.c
  param/param.c
  collect/collect.c
  report/report_help.c
  report/stat_report.c
  report/json_report.c
  report/mqc_report.c
  report/outfile.c
  ${CMAKE_CURRENT_BINARY_DIR}/convert_tables.h
  )

if(Threads_FOUND)
  list(APPEND source_files thread/thread_data.c)
endif(Threads_FOUND)

add_executable(samstat)

target_sources(samstat
  PRIVATE
  ${source_files}
  )

target_include_directories(samstat PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

if(SAMSTAT_USE_SYSTEM_HTSLIB)
  target_link_libraries(samstat PRIVATE PkgConfig::HTSLIB)
else()
  add_dependencies(samstat htslib)
  target_include_directories(samstat PRIVATE ${HTSLIB_SRC_DIR})
  target_link_directories(samstat PRIVATE  ${HTSLIB_SRC_DIR})
  target_link_libraries(samstat PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)
endif()

if(Threads_FOUND)
  target_link_libraries(samstat PRIVATE Threads::Threads)
endif()

target_link_libraries(samstat PRIVATE tldlite module_plot)

if(SAMSTAT_SANITIZE)
  target_compile_options(samstat PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(samstat PRIVATE -fsanitize=address,undefined)
endif()
