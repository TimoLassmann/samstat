add_subdirectory(util)
add_subdirectory(plot)

# collect.c includes collect/convert_tables.h, which is committed to the source tree: it used to
# be written at build time by converttablegen (convert_tables.c), which cannot run on the build
# host of a cross build. The test build:convert_tables_in_sync builds and runs the generator and
# fails if its output differs from the committed header.

set(source_files
  samstat.c
  htsinterface/htsglue.c
  sambamparse/sam_bam_parse.c
  param/param.c
  collect/collect.c
  report/report_help.c
  report/stat_report.c
  report/json_report.c
  report/mqc_report.c
  thread/thread_data.c
  )

add_executable(samstat)

target_sources(samstat
  PRIVATE
  ${source_files}
  )

# htslib needs pthreads in any case (htslib.cmake requires them), so samstat always runs its
# reader and workers in threads; there is no build without them.
target_compile_definitions(samstat PRIVATE HAVE_PTHREADS)

target_link_libraries(samstat PRIVATE htslib::hts Threads::Threads tldlite module_plot)

if(SAMSTAT_SANITIZE)
  target_compile_options(samstat PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(samstat PRIVATE -fsanitize=address,undefined)
endif()
