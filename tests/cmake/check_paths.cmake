# Fails if any path below ROOT (<build>/tests/out) is longer than MAX_LEN characters or more
# than MAX_DEPTH levels deep, counted from ROOT. A deep or very long path left in the build
# tree breaks backup tools (Time Machine) and macOS's own temp-folder cleanup, so this test
# runs after all others and names each offending path (tests/README.md, rule 5).
#
#   cmake -DROOT=<dir> -DMAX_LEN=<n> -DMAX_DEPTH=<n> -P check_paths.cmake

if(NOT IS_DIRECTORY "${ROOT}")
  return()
endif()

file(GLOB_RECURSE paths LIST_DIRECTORIES true RELATIVE "${ROOT}" "${ROOT}/*")
set(bad "")
foreach(p IN LISTS paths)
  string(LENGTH "${p}" len)
  string(REGEX MATCHALL "/" slashes "${p}")
  list(LENGTH slashes depth)
  math(EXPR depth "${depth} + 1")
  if(len GREATER MAX_LEN OR depth GREATER MAX_DEPTH)
    string(APPEND bad "  ${len} characters, ${depth} levels: ${p}\n")
  endif()
endforeach()

if(bad)
  message(FATAL_ERROR "paths below ${ROOT} longer than ${MAX_LEN} characters or deeper "
                      "than ${MAX_DEPTH} levels:\n${bad}")
endif()
