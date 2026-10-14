# Compare a samstat JSON report against a golden file, ignoring the 'meta' member
# (version and file name are the only things allowed to differ).
# Usage: cmake -DACTUAL=<json> -DEXPECTED=<json> -P compare_golden.cmake
#        cmake -DACTUAL=<json> -DEXPECTED=<json> -DUPDATE=ON -P compare_golden.cmake
# The second form overwrites EXPECTED with ACTUAL (regenerate goldens).

cmake_minimum_required(VERSION 3.19)

if(NOT ACTUAL OR NOT EXISTS "${ACTUAL}")
  message(FATAL_ERROR "compare_golden: ACTUAL not found: '${ACTUAL}'")
endif()
if(NOT EXPECTED)
  message(FATAL_ERROR "compare_golden: EXPECTED not given")
endif()

if(UPDATE)
  get_filename_component(dir "${EXPECTED}" DIRECTORY)
  file(MAKE_DIRECTORY "${dir}")
  file(COPY_FILE "${ACTUAL}" "${EXPECTED}")
  message(STATUS "compare_golden: updated ${EXPECTED}")
  return()
endif()

if(NOT EXISTS "${EXPECTED}")
  message(FATAL_ERROR "compare_golden: EXPECTED not found: ${EXPECTED}\n(create it with -DUPDATE=ON)")
endif()

file(READ "${ACTUAL}" actual)
file(READ "${EXPECTED}" expected)

string(JSON actual ERROR_VARIABLE err REMOVE "${actual}" meta)
if(NOT err STREQUAL "NOTFOUND")
  message(FATAL_ERROR "compare_golden: ACTUAL is not valid JSON or has no 'meta': ${err}")
endif()
string(JSON expected ERROR_VARIABLE err REMOVE "${expected}" meta)
if(NOT err STREQUAL "NOTFOUND")
  message(FATAL_ERROR "compare_golden: EXPECTED is not valid JSON or has no 'meta': ${err}")
endif()

string(JSON same EQUAL "${actual}" "${expected}")
if(same)
  return()
endif()

# Report the first differing top-level key.
string(JSON n_act LENGTH "${actual}")
string(JSON n_exp LENGTH "${expected}")
set(first "")
math(EXPR last "${n_act}-1")
foreach(i RANGE ${last})
  string(JSON key MEMBER "${actual}" ${i})
  string(JSON e ERROR_VARIABLE err GET "${expected}" "${key}")
  if(NOT err STREQUAL "NOTFOUND")
    set(first "${key} (missing in EXPECTED)")
    break()
  endif()
  # compare the member values as JSON documents (all non-meta members are numbers, arrays or objects)
  string(JSON ja GET "${actual}" "${key}")
  string(JSON je GET "${expected}" "${key}")
  string(JSON member_same EQUAL "${ja}" "${je}")
  if(NOT member_same)
    set(first "${key}")
    break()
  endif()
endforeach()
if(first STREQUAL "")
  set(first "(a key present only in EXPECTED; ${n_act} vs ${n_exp} members)")
endif()
message(FATAL_ERROR "compare_golden: ${ACTUAL} differs from ${EXPECTED}\nfirst differing top-level key: ${first}\n(regenerate with -DUPDATE=ON if the change is intended)")
