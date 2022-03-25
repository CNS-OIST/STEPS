add_library(stepsrng STATIC
    rng.cpp
    mt19937.cpp
    r123.cpp
    create.cpp
)

set_property(TARGET stepsrng PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsrng PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "9.0"
   OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
  set_source_files_properties(r123.cpp
                              create.cpp
                              PROPERTIES
                              COMPILE_FLAGS
                              "-Wno-old-style-cast -Wno-expansion-to-defined")
else()
  set_source_files_properties(r123.cpp
                              create.cpp
                              PROPERTIES
                              COMPILE_FLAGS
                              "-Wno-old-style-cast")
endif()

target_link_libraries(stepsrng stepsutil stepsmath)
