# tld is a submodule of the upstream repository (one tld codebase, no vendored copy).
# Required: -DSOURCE_DIR=<repo>
file(READ ${SOURCE_DIR}/.gitmodules mods)
string(TOLOWER "${mods}" mods)
if(NOT mods MATCHES "path = src/tld" OR NOT mods MATCHES "url = https://github\\.com/timolassmann/tld(\\.git)?")
  message(FATAL_ERROR ".gitmodules does not point src/tld at github.com/timolassmann/tld:\n${mods}")
endif()
if(NOT EXISTS ${SOURCE_DIR}/src/tld/CMakeLists.txt)
  message(FATAL_ERROR "src/tld is empty: git submodule update --init --recursive")
endif()
find_program(GIT git)
if(GIT AND EXISTS ${SOURCE_DIR}/.git)
  execute_process(COMMAND ${GIT} -C ${SOURCE_DIR} ls-tree HEAD src/tld
    OUTPUT_VARIABLE tree RESULT_VARIABLE rc OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(rc EQUAL 0 AND NOT tree MATCHES "^160000 commit [0-9a-f]+\tsrc/tld$")
    message(FATAL_ERROR "src/tld is not a pinned submodule commit: '${tree}'")
  endif()
endif()
