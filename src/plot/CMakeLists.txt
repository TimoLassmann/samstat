# SVG primitives for the report (WP1.6). A separate library so that the
# sanitizer build instruments it and the unit test can link it on its own.
add_library(svg_plot STATIC svg_plot.c)
target_include_directories(svg_plot INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(svg_plot PRIVATE tld-dev)
if(UNIX)
  target_link_libraries(svg_plot PRIVATE m)
endif()
if(SAMSTAT_SANITIZE)
  target_compile_options(svg_plot PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
endif()

add_executable(svg_plot_test svg_plot_test.c)
target_link_libraries(svg_plot_test PRIVATE tld-dev svg_plot)
if(UNIX)
  target_link_libraries(svg_plot_test PRIVATE m)
endif()
if(SAMSTAT_SANITIZE)
  target_compile_options(svg_plot_test PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(svg_plot_test PRIVATE -fsanitize=address,undefined)
endif()

# report_plots.c (WP4.2b) draws the six standard plots from a finalised stat_collection. It uses
# the wording of report/report_help.c and the tables of collect/collect.c, which the final
# executable links in; module_plot itself does not need them until report_plots.o is pulled.
add_library(module_plot STATIC
  plot.c
  plot_interval.c
  plot_group.c
  report_plots.c
  )

target_include_directories(module_plot
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_plot PUBLIC svg_plot)
target_link_libraries(module_plot PRIVATE tld-dev)

add_executable(unit-plot plot_test.c)
target_link_libraries(unit-plot PRIVATE tld-dev module_plot)

# add_test(
#     NAME unit-str
#     COMMAND unit-str
# )
