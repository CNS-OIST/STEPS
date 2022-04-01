# -----------------------------------------------------------------------------
# * Find Sundials includes and libraries.
#
# This module finds if Sundials is installed and determines where the include
# files and libraries are.  This code sets the following variables:
#
# * SUNDIALS_FOUND         = Sundials was found
# * SUNDIALS_LIBRARY_DIR   = path to where libraries can be found
# * SUNDIALS_INCLUDE_DIR   = path to where header files can be found
# * SUNDIALS_LIBRARIES     = link libraries for Sundials
# -----------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)

find_path(SUNDIALS_DIR include/sundials/sundials_config.h
          HINTS ${SUNDIALS_INSTALL_DIR} ENV
          PATHS $ENV{HOME}/sundials
          DOC "Sundials Directory")

message(STATUS "SUNDIALS_DIR=${SUNDIALS_DIR}")

if(SUNDIALS_DIR)
  set(SUNDIALS_FOUND YES)

  set(SUNDIALS_INCLUDE_DIR ${SUNDIALS_DIR}/include)
  set(SUNDIALS_LIBRARY_DIR ${SUNDIALS_DIR}/lib)

  # The set of required Sundials libraries
  set(SUNDIALS_REQUIRED_LIBS sundials_nvecserial sundials_ida sundials_kinsol)
  if(OPENMP_ENABLE)
    list(APPEND SUNDIALS_REQUIRED_LIBS sundials_nvecopenmp)
  endif(OPENMP_ENABLE)

  if(LOAD_CVODE)
    list(APPEND SUNDIALS_REQUIRED_LIBS sundials_cvode)
  endif(LOAD_CVODE)

  if(LOAD_ARKODE)
    list(APPEND SUNDIALS_REQUIRED_LIBS sundials_arkode)
  endif(LOAD_ARKODE)

  foreach(TEST_LIB IN LISTS SUNDIALS_REQUIRED_LIBS)

    message(STATUS "Checking for TEST_LIB=${TEST_LIB}")

    # Need to make sure variable to search for isn't set
    unset(SUNDIALS_LIB CACHE)

    find_library(SUNDIALS_LIB
                 NAMES ${TEST_LIB}
                 HINTS ${SUNDIALS_LIBRARY_DIR})

    message(STATUS "SUNDIALS_LIB=${SUNDIALS_LIB}")
    if(SUNDIALS_LIB)
      list(APPEND SUNDIALS_LIBRARIES ${SUNDIALS_LIB})
    else(SUNDIALS_LIB)
      message(
        FATAL_ERROR "Could not find required Sundials library : ${TEST_LIB}")
    endif(SUNDIALS_LIB)

  endforeach(TEST_LIB)

  message(STATUS "CMAKE_FIND_LIBRARY_SUFFIXES=${CMAKE_FIND_LIBRARY_SUFFIXES}")
  message(STATUS "SUNDIALS_FOUND=${SUNDIALS_FOUND}")
  message(STATUS "SUNDIALS_LIBRARY_DIR=${SUNDIALS_LIBRARY_DIR}")
  message(STATUS "SUNDIALS_INCLUDE_DIR=${SUNDIALS_INCLUDE_DIR}")
  message(STATUS "SUNDIALS_LIBRARIES=${SUNDIALS_LIBRARIES}")

  # Extract version
  set(SUNDIALS_GET_VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/sundials_get_version.cpp")
  file(WRITE ${SUNDIALS_GET_VERSION_FILE}
    "
    #include <iostream>
    #include <sundials/sundials_config.h>
    #ifndef SUNDIALS_VERSION
    # define SUNDIALS_VERSION SUNDIALS_PACKAGE_VERSION
    #endif
    int main() {
        std::cout << SUNDIALS_VERSION;
    }")

  try_run(SUNDIALS_RUN_RESULT
    SUNDIALS_COMPILE_RESULT
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SUNDIALS_GET_VERSION_FILE}
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${SUNDIALS_INCLUDE_DIR}"
    COMPILE_OUTPUT_VARIABLE SUNDIALS_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE SUNDIALS_VERSION
    )

  if(NOT SUNDIALS_COMPILE_RESULT)
    message(SEND_ERROR "Could not compile ${SUNDIALS_GET_VERSION_FILE}:\n${SUNDIALS_COMPILE_OUTPUT}")
  else()
    # Cleanup
    file(REMOVE ${SUNDIALS_GET_VERSION_FILE})

    message(STATUS "SUNDIALS_VERSION=${SUNDIALS_VERSION}")

    # Split version
    string(REPLACE "." ";" VERSION_LIST ${SUNDIALS_VERSION})
    list(GET VERSION_LIST 0 SUNDIALS_VERSION_MAJOR)
    list(GET VERSION_LIST 1 SUNDIALS_VERSION_MINOR)
    list(GET VERSION_LIST 2 SUNDIALS_VERSION_PATCH)
  endif()

  if(CMAKE_VERSION VERSION_LESS 3.12.0)
    add_definitions(-DSTEPS_SUNDIALS_VERSION_MAJOR=${SUNDIALS_VERSION_MAJOR})
  else()
    add_compile_definitions(STEPS_SUNDIALS_VERSION_MAJOR=${SUNDIALS_VERSION_MAJOR})
  endif()

else(SUNDIALS_DIR)
  set(SUNDIALS_FOUND NO)
endif(SUNDIALS_DIR)

if(CMAKE_VERSION VERSION_LESS 3.19.0)
  find_package_handle_standard_args(SUNDIALS
                                    FOUND_VAR SUNDIALS_FOUND
                                    REQUIRED_VARS SUNDIALS_LIBRARIES SUNDIALS_INCLUDE_DIR
                                    VERSION_VAR SUNDIALS_VERSION)
else()
  find_package_handle_standard_args(SUNDIALS
                                    FOUND_VAR SUNDIALS_FOUND
                                    REQUIRED_VARS SUNDIALS_LIBRARIES SUNDIALS_INCLUDE_DIR
                                    VERSION_VAR SUNDIALS_VERSION
                                    HANDLE_VERSION_RANGE)
endif()