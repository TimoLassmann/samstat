add_library(module_plot STATIC
  plot.c
  plot_interval.c
  plot_group.c
  )

target_include_directories(module_plot 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_plot PRIVATE tld-dev)

add_executable(unit-plot plot_test.c) 
target_link_libraries(unit-plot PRIVATE tld-dev module_plot)

# add_test(
#     NAME unit-str
#     COMMAND unit-str
# )
