# Check values inside a JSON file.
# Usage: cmake -DJSON=<file> -DCHECKS="path1=value1;path2=value2" -P check_json.cmake
#
# A path is a dot separated list of object keys and array indices, e.g.
#   summary.total_reads       read1.length.mapq_ge30.42
# 'path#len=N' checks the number of elements of an array (or members of an object).
# 'path=value' compares the scalar at path; booleans are written true / false.
# 'path>=N' and 'path<=N' compare a numeric scalar (or, with #len, a length).
# With no operator the check only requires that the path exists.
# Values must not contain ';' or '='. Path keys must not contain '<' or '>' before '='.

cmake_minimum_required(VERSION 3.19)

if(NOT JSON OR NOT EXISTS "${JSON}")
  message(FATAL_ERROR "check_json: JSON file not found: '${JSON}'")
endif()
if(NOT DEFINED CHECKS)
  message(FATAL_ERROR "check_json: no CHECKS given")
endif()

# add_test escapes ";" as backslash-semicolon - undo that so CHECKS is a proper list
string(REPLACE "\\;" ";" CHECKS "${CHECKS}")

file(READ "${JSON}" doc)

# Turn 'a.b.0.c' into the arguments of string(JSON GET): 'a;b;0;c'
function(json_path_to_list path out)
  string(REPLACE "." ";" parts "${path}")
  set(${out} "${parts}" PARENT_SCOPE)
endfunction()

set(n_fail 0)
set(n_done 0)
foreach(check IN LISTS CHECKS)
  if(check STREQUAL "")
    continue()
  endif()
  math(EXPR n_done "${n_done}+1")

  # operator: first of '>=', '<=', '=' (the two-character operators contain '=')
  set(op "")
  string(FIND "${check}" "=" eq)
  if(NOT eq EQUAL -1)
    set(op "=")
    set(oplen 1)
    set(opat ${eq})
    if(eq GREATER 0)
      math(EXPR prev "${eq}-1")
      string(SUBSTRING "${check}" ${prev} 1 prevchar)
      if(prevchar STREQUAL ">" OR prevchar STREQUAL "<")
        set(op "${prevchar}=")
        set(oplen 2)
        set(opat ${prev})
      endif()
    endif()
  endif()
  if(op STREQUAL "")
    set(path "${check}")
    set(has_expected FALSE)
    set(expected "")
  else()
    string(SUBSTRING "${check}" 0 ${opat} path)
    math(EXPR vstart "${opat}+${oplen}")
    string(SUBSTRING "${check}" ${vstart} -1 expected)
    set(has_expected TRUE)
  endif()

  set(want_len FALSE)
  if(path MATCHES "^(.*)#len$")
    set(path "${CMAKE_MATCH_1}")
    set(want_len TRUE)
  endif()

  json_path_to_list("${path}" keys)
  if(want_len)
    string(JSON actual ERROR_VARIABLE err LENGTH "${doc}" ${keys})
  else()
    string(JSON actual ERROR_VARIABLE err GET "${doc}" ${keys})
  endif()

  if(NOT err STREQUAL "NOTFOUND")
    message(SEND_ERROR "check_json: '${check}': cannot resolve path '${path}' in ${JSON}: ${err}")
    math(EXPR n_fail "${n_fail}+1")
    continue()
  endif()

  if(has_expected AND (op STREQUAL ">=" OR op STREQUAL "<="))
    if(NOT actual MATCHES "^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$" OR NOT expected MATCHES "^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$")
      message(SEND_ERROR "check_json: '${check}': non-numeric operand (got '${actual}', limit '${expected}') (${JSON})")
      math(EXPR n_fail "${n_fail}+1")
    elseif(op STREQUAL ">=" AND actual LESS expected)
      message(SEND_ERROR "check_json: mismatch for '${path}': expected >= ${expected}, got '${actual}' (${JSON})")
      math(EXPR n_fail "${n_fail}+1")
    elseif(op STREQUAL "<=" AND actual GREATER expected)
      message(SEND_ERROR "check_json: mismatch for '${path}': expected <= ${expected}, got '${actual}' (${JSON})")
      math(EXPR n_fail "${n_fail}+1")
    endif()
  elseif(has_expected)
    # string(JSON GET) reports booleans as ON / OFF
    if(expected STREQUAL "true")
      set(expected "ON")
    elseif(expected STREQUAL "false")
      set(expected "OFF")
    endif()
    if(NOT actual STREQUAL expected)
      message(SEND_ERROR "check_json: mismatch for '${path}': expected '${expected}', got '${actual}' (${JSON})")
      math(EXPR n_fail "${n_fail}+1")
    endif()
  endif()
endforeach()

if(n_done EQUAL 0)
  message(FATAL_ERROR "check_json: CHECKS contains no check (a typo in the test?)")
endif()
if(n_fail GREATER 0)
  message(FATAL_ERROR "check_json: ${n_fail} of ${n_done} checks failed")
endif()
message(STATUS "check_json: ${n_done} checks passed")
