# tldlite: the parts of the tld library samstat uses, as one header (tldlite.h) and the
# translation unit that compiles it. Every target that includes tldlite.h links this.
add_library(tldlite STATIC tldlite.c)
target_include_directories(tldlite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# the math library, for every target that includes tldlite.h (which includes math.h)
if(UNIX)
  target_link_libraries(tldlite PUBLIC m)
endif()
if(SAMSTAT_SANITIZE)
  target_compile_options(tldlite PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  # every executable that links tldlite then links the sanitizer runtime too
  target_link_options(tldlite PUBLIC -fsanitize=address,undefined)
endif()
