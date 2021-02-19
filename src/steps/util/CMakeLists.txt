add_library(
    stepsutil
    checkid.cpp
    ../init.cpp
    ../finish.cpp
    ../error.cpp
    ${EASYLOGGINGPP_INCLUDE_DIRS}/easylogging++.cc
)

set_property(TARGET stepsutil PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsutil
    PUBLIC .
    PUBLIC ../../)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(
    easyloggingpp_flags
    "-Wno-double-promotion -Wno-shadow -Wno-unused-parameter -Wno-useless-cast \
    -Wno-null-dereference"
    )
else()
  set(
    easyloggingpp_flags
    "-Wno-double-promotion -Wno-shadow -Wno-unused-parameter -Wno-null-dereference \
    -Wno-range-loop-analysis"
    )
endif()
set_source_files_properties("${EASYLOGGINGPP_INCLUDE_DIRS}/easylogging++.cc"
                            PROPERTIES
                            COMPILE_FLAGS
                            ${easyloggingpp_flags})