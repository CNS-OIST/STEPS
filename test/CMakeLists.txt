set(gtest_src "${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest/googletest")
set(gtest_bin "${CMAKE_CURRENT_BINARY_DIR}/third_party/gtest")
set(catch2 ${CMAKE_CURRENT_SOURCE_DIR}/third_party/catch2)

if(NOT EXISTS ${gtest_src}/CMakeLists.txt)
  message(
    FATAL_ERROR
      " Submodules do not seems to be initialized \n Please run 'git submodule update --init --recursive'"
    )
endif()

add_subdirectory(${gtest_src} ${gtest_bin} EXCLUDE_FROM_ALL)

set(inc_dirs "${CMAKE_SOURCE_DIR}/src" "${gtest_src}/include" "${catch2}/include" )
set(libs libsteps gtest_main)

include_directories(${inc_dirs})

# Check Python dependencies
set(PYMODS matplotlib numpy scipy)
set(PYMODS_FOUND "TRUE")
include(FindPythonModule)
foreach(MOD ${PYMODS})
  find_python_module(${MOD})
  string(TOUPPER ${MOD} MODU)
  if(NOT ${MODU}_FOUND)
    set(PYMODS_FOUND "FALSE")
  endif()
endforeach()
if(NOT PYMODS_FOUND)
  message(WARNING "One or more python module is missing, some tests might fail.")
  message(STATUS "Consider installing the missing python module or disabling tests \"cmake -DBUILD_TESTING=OFF ..\".")
endif()

add_subdirectory(unit)
add_subdirectory(validation)
