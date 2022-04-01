if(NOT EXISTS "${EASYLOGGINGPP_INCLUDE_DIRS}/easylogging++.cc")
  message(
    FATAL_ERROR
      "Submodule easyloggingpp does not seem to be initialized \n Please run 'git submodule update --init --recursive'"
    )
endif()


add_subdirectory(tracker)
add_subdirectory(profile)

add_library(stepsutil STATIC
    checkid.cpp
    debug.cpp
    error.cpp
    finish.cpp
    init.cpp
    ${EASYLOGGINGPP_INCLUDE_DIRS}/easylogging++.cc
)

set_property(TARGET stepsutil PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsutil PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(
    easyloggingpp_flags
    "-Wno-double-promotion -Wno-shadow -Wno-unused-parameter -Wno-useless-cast \
    -Wno-null-dereference")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "11")
    set(easyloggingpp_flags "${easyloggingpp_flags} -Wno-range-loop-construct")
  endif()
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

add_definitions(-DELPP_NO_DEFAULT_LOG_FILE=1 -DELPP_STL_LOGGING=1
                -DELPP_DISABLE_DEFAULT_CRASH_HANDLING=1)

target_link_libraries(stepsutil PUBLIC stepstracker instrumentor)
if(PETSC_FOUND AND USE_PETSC)
  target_link_libraries(stepsutil PUBLIC PETSC::PETSC)
endif()

if(STEPS_USE_DIST_MESH)
    target_link_libraries(stepsutil PUBLIC Omega_h::omega_h "${Gmsh_LIBRARIES}")
endif()
