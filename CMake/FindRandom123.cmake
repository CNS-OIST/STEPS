#
# Find Random123 includes
#
# Random123 It can be found at:
#
# * Random123_INCLUDE_DIR - where to find Random123.h
# * Random123_FOUND       - boolean indicating if Random123 was found.

find_path(Random123_INCLUDE_DIR Random123/threefry.h Random123/MicroURNG.hpp)

if(Random123_INCLUDE_DIR)
  set(Random123_FOUND "YES")
endif(Random123_INCLUDE_DIR)

if(Random123_FOUND)
  if(NOT Random123_FIND_QUIETLY)
    message(STATUS "Found Random123:${Random123_INCLUDE_DIR}")
  endif(NOT Random123_FIND_QUIETLY)
else(Random123_FOUND)
  if(Random123_FIND_REQUIRED)
    message(FATAL_ERROR "Random123 not found!")
  endif(Random123_FIND_REQUIRED)
endif(Random123_FOUND)
