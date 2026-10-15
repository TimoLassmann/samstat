# tldlite.h, the parts of the tld library samstat uses (label: util).
add_executable(unit_tldlite unit_tldlite.c)
target_link_libraries(unit_tldlite PRIVATE tldlite)
if(SAMSTAT_SANITIZE)
  target_compile_options(unit_tldlite PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(unit_tldlite PRIVATE -fsanitize=address,undefined)
endif()
add_test(NAME util:unit_tldlite COMMAND unit_tldlite)
set_tests_properties(util:unit_tldlite PROPERTIES LABELS util)
