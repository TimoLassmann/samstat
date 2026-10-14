# SVG primitive tests (label: report).
#
# svg:unit runs src/plot/svg_plot_test: ticks and labels, every primitive,
# byte-for-byte determinism, and a well-formedness check of every SVG produced
# (tag balance, attribute quoting, role/title/desc, no non-finite numbers).
# The test renders into memory only; nothing is written to disk.

add_test(NAME svg:unit COMMAND svg_plot_test)
set_tests_properties(svg:unit PROPERTIES LABELS report)

# svg:report_plots runs src/plot/report_plots_test (WP4.2b): every plot renders for synthetic
# finalised collections (short, long, paired, unpaired, unaligned, no qualities, no MD, nothing
# mapped, empty), the markup is well formed, Read 1 and Read 2 share their axes, the bytes are
# deterministic and locale independent, the worst case (paired 999 bp) stays small, and the base
# and MAPQ palettes share no hex value. The test links the collector (for the length tables and the
# stat_collection layout) the same way collect:unit does, and the report wording (report_help.c).
set(SVG_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${SVG_GEN_DIR})
add_custom_command(
  OUTPUT ${SVG_GEN_DIR}/convert_tables.h
  COMMAND converttablegen ${SVG_GEN_DIR}/convert_tables.h
  DEPENDS converttablegen)
add_executable(report_plots_test
  ${PROJECT_SOURCE_DIR}/src/plot/report_plots_test.c
  ${PROJECT_SOURCE_DIR}/src/report/report_help.c
  ${PROJECT_SOURCE_DIR}/src/collect/collect.c
  ${PROJECT_SOURCE_DIR}/src/htsinterface/htsglue.c
  ${PROJECT_SOURCE_DIR}/src/sambamparse/sam_bam_parse.c
  ${SVG_GEN_DIR}/convert_tables.h)
target_include_directories(report_plots_test PRIVATE ${SVG_GEN_DIR})
if(SAMSTAT_USE_SYSTEM_HTSLIB)
  target_link_libraries(report_plots_test PRIVATE PkgConfig::HTSLIB)
else()
  add_dependencies(report_plots_test htslib)
  target_include_directories(report_plots_test PRIVATE ${HTSLIB_SRC_DIR})
  target_link_libraries(report_plots_test PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)
endif()
if(Threads_FOUND)
  target_link_libraries(report_plots_test PRIVATE Threads::Threads)
endif()
target_link_libraries(report_plots_test PRIVATE module_plot tld-dev)
if(UNIX)
  target_link_libraries(report_plots_test PRIVATE m)
endif()
if(SAMSTAT_SANITIZE)
  target_compile_options(report_plots_test PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(report_plots_test PRIVATE -fsanitize=address,undefined)
endif()

add_test(NAME svg:report_plots COMMAND report_plots_test)
set_tests_properties(svg:report_plots PROPERTIES LABELS report)
