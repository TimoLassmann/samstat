# SVG primitive tests (label: report).
#
# svg:unit runs src/plot/svg_plot_test: ticks and labels, every primitive,
# byte-for-byte determinism, and a well-formedness check of every SVG produced
# (tag balance, attribute quoting, role/title/desc, no non-finite numbers).
# The test renders into memory only; nothing is written to disk.

add_test(NAME svg:unit COMMAND svg_plot_test)
set_tests_properties(svg:unit PROPERTIES LABELS report)
