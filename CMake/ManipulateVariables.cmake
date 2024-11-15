# Routines for finding variables and adding or removing words in their values.

function(match_variables return regex)
  set(_all_vars)
  set(_matched_vars)
  get_cmake_property(_all_vars VARIABLES)
  foreach(_var ${_all_vars})
    set(_match)
    string(REGEX MATCH "${regex}" _match ${_var})
    if(_match)
      list(APPEND _matched_vars ${_var})
    endif()
  endforeach()

  set("${return}"
      "${_matched_vars}"
      PARENT_SCOPE)
endfunction()

function(string_split return input)
  if("${ARGV2}" STREQUAL "")
    set(_sep " +")
  else()
    set(_sep ${ARGV2})
  endif()

  string(REGEX REPLACE "${_sep}" ";" _out "${input}")
  set("${return}"
      "${_out}"
      PARENT_SCOPE)
endfunction()

function(string_join return separator)
  string(REGEX REPLACE ";" "${separator}" _out "${ARGN}")
  set("${return}"
      "${_out}"
      PARENT_SCOPE)
endfunction()

function(remove_word var word)
  if(NOT "${${var}}" STREQUAL "")
    string_split(_words "${${var}}")
    list(REMOVE_ITEM _words "${word}")
    string_join(_new_value " " ${_words})
    set("${var}"
        ${_new_value}
        PARENT_SCOPE)
  endif()
endfunction()
