# ~~~
# For older CMake ( < 3.18 ) libraries found through FindBLA do not
# define a library target that carries all dependencies. This function
# supplies the missing target by parsing the relevant lists.
#
# current implementation processes
# - shared libs
# - linking flags (-l)
# ~~~

function(add_library_target)
  set(oneValueArgs NAME)
  set(multiValueArgs LIBRARIES INCLUDE_DIRECTORIES)
  set(options NAME_LOWER_CASE)
  cmake_parse_arguments(ADT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(_base_name "${ADT_NAME}::")
  if(ADT_NAME_LOWER_CASE)
    string(TOLOWER "${ADT_NAME}" lw_name)
    set(_lib_name "${_base_name}${lw_name}")
  else()
    set(_lib_name "${_base_name}${ADT_NAME}")
  endif()

  # check if modern target exists
  if(TARGET "${_lib_name}")
    message(STATUS "Modern target ${_lib_name} already exists, skipping setup")
  else()
    message(STATUS "Added modern target ${_lib_name}")
    add_library("${_lib_name}" INTERFACE IMPORTED)
    # some FindBLA are buggy...
    list(REMOVE_DUPLICATES ADT_LIBRARIES)
    foreach(_chunk ${ADT_LIBRARIES})
      # process shared libs
      set(_is_shared ${_chunk})
      list(FILTER _is_shared INCLUDE REGEX "\\${CMAKE_SHARED_LIBRARY_SUFFIX}")
      if(_is_shared)
        get_filename_component(_lib ${_chunk} NAME)
        get_filename_component(_name ${_chunk} NAME_WE)
        if(NOT TARGET "${_base_name}${_name}")
          add_library("${_base_name}${_name}" SHARED IMPORTED GLOBAL)
          set_property(TARGET "${_base_name}${_name}" PROPERTY IMPORTED_LOCATION "${_chunk}")
          list(APPEND _all_so_libs "${_base_name}${_name}")
        endif()
      endif()
      # process -l prepended libs
      string(REGEX MATCH "(^-l)(.*)" _tmp ${_chunk})
      if(CMAKE_MATCH_0)
        list(APPEND _all_raw_libs "${CMAKE_MATCH_2}")
      endif()
      # process libs with no extension and no -l
      if(NOT CMAKE_MATCH_0 AND NOT _is_shared)
        list(APPEND _all_raw_libs "${_chunk}")
      endif()
    endforeach()
    set_property(TARGET "${_lib_name}" PROPERTY INTERFACE_LINK_LIBRARIES
                                                "${_all_so_libs};${_all_raw_libs}")
    set_property(TARGET "${_lib_name}" PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                "${ADT_INCLUDE_DIRECTORIES}")
  endif()
endfunction()
