add_library(
  stepsrng STATIC
  create.cpp
  create.hpp
  mt19937.cpp
  mt19937.hpp
  r123.cpp
  r123.hpp
  rng.cpp
  rng.hpp
  std_mt19937.cpp
  std_mt19937.hpp)

set_property(TARGET stepsrng PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsrng PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "9.0"
   OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
  set_source_files_properties(
    r123.cpp create.cpp PROPERTIES COMPILE_FLAGS "-Wno-old-style-cast -Wno-expansion-to-defined")
else()
  set_source_files_properties(r123.cpp create.cpp PROPERTIES COMPILE_FLAGS "-Wno-old-style-cast")
endif()

target_link_libraries(stepsrng stepsutil stepsmath)
