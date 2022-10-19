add_library(plot STATIC
  plot.c
  plot_interval.c
  )

target_include_directories(plot 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(plot PRIVATE tld-dev)

add_executable(unit-plot plot_test.c) 
target_link_libraries(unit-plot PRIVATE tld-dev plot)

# add_test(
#     NAME unit-str
#     COMMAND unit-str
# )
