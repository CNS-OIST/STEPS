cmake_minimum_required(VERSION 3.10.2)
if(POLICY CMP0042)
  # Use rpath on Mac OS X
  cmake_policy(SET CMP0042 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

project(STEPS)

set(VERSION_MAJOR 3)
set(VERSION_MINOR 6)
set(VERSION_PATCH 0)
set(STEPS_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# OPTIONS
option(BUILD_STOCHASTIC_TESTS "Build stochastic tests" ON)
option(TARGET_NATIVE_ARCH "Generate non-portable arch-specific code" ON)
option(USE_BDSYSTEM_LAPACK "Use new BDSystem/Lapack code for E-Field solver"
       OFF)

option(USE_BUNDLE_EASYLOGGINGPP "Use bundled version of easylogging" ON)
option(USE_BUNDLE_RANDOM123 "Use bundled version of random123" ON)
option(USE_BUNDLE_SUNDIALS "Use bundled version of cvode" ON)
option(USE_64_BITS_INDICES "Use 64bits indices instead of 32" OFF)

set(USE_MPI "Default" CACHE STRING "Use MPI for parallel solvers")
# Toggleable options
set_property(CACHE USE_MPI PROPERTY STRINGS "Default;True;False")
set(USE_PETSC "Default"
    CACHE STRING "Use PETSC library for parallel E-Field solver")
set_property(CACHE USE_PETSC PROPERTY STRINGS "Default;True;False")

if(USE_BUNDLE_EASYLOGGINGPP)
  set(EASYLOGGINGPP_INCLUDE_DIRS
      "${CMAKE_SOURCE_DIR}/src/third_party/easyloggingpp/src")
else()
  find_package(PkgConfig)
  pkg_check_modules(EASYLOGGINGPP REQUIRED easyloggingpp)
endif()
include_directories(SYSTEM ${EASYLOGGINGPP_INCLUDE_DIRS})

if(USE_BUNDLE_RANDOM123)
  include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/src/third_party")
else()
  find_package(Random123 REQUIRED)
  include_directories(SYSTEM ${Random123_INCLUDE_DIR})
endif()

if(USE_64_BITS_INDICES)
  add_definitions(-DSTEPS_USE_64BITS_INDICES=1)
endif()

if(USE_BUNDLE_SUNDIALS)
  set(cvode "${CMAKE_SOURCE_DIR}/src/third_party/cvode-2.6.0/src")
  set(cvode_include "${CMAKE_SOURCE_DIR}/src/third_party/cvode-2.6.0/include")
  include_directories(SYSTEM ${cvode} ${cvode_include})
else()
  set(LOAD_CVODE ON)
  find_package(SUNDIALS REQUIRED)
  include_directories(SYSTEM ${SUNDIALS_INCLUDE_DIR})
  link_directories(${SUNDIALS_LIBRARY_DIR})
endif()

# Valgrind
set(VALGRIND "" CACHE STRING "Valgrind plus arguments for testing")
if(NOT VALGRIND STREQUAL "")
  file(GLOB valgrind_supps ${CMAKE_CURRENT_SOURCE_DIR}/test/ci/valgrind/*.supp)
  foreach(valgrind_supp IN LISTS valgrind_supps)
    list(APPEND valgrind_supps_cmd --suppressions=${valgrind_supp})
  endforeach()
  set(VALGRIND "${VALGRIND}" ${valgrind_supps_cmd})
endif()

# Compiler options
# -------------------------------------------------------------------
include(CMake/Compiler.cmake)

include(CMake/steps_portability_check.cmake)

option(USE_CLANG_TIDY "Perform C++ static analysis while compiling" FALSE)
if(USE_CLANG_TIDY)
  find_package(CLANG_TIDY REQUIRED)
  if(CLANG_TIDY_FOUND)
    cmake_minimum_required(VERSION 3.6)
    set(CLANG_TIDY_ARGS "" CACHE STRING "clang-tidy command options")
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}" ${CLANG_TIDY_ARGS})
  else()
    message(ERROR "Could not find clang-tidy utility")
  endif()
endif()

# Nuke any -DNDEBUG in the compiler options introduced by CMake.
include(CMake/ManipulateVariables.cmake)
foreach(var_name ${flag_vars})
  remove_word(${var_name} "-DNDEBUG")
endforeach()

# for debug symbol
include(CMake/UseDebugSymbols.cmake)

add_options(ALL_LANGUAGES
            ALL_BUILDS
            "-Wall"
            "-Wextra"
            "-Wshadow"
            "-Wnon-virtual-dtor"
            "-Wunused"
            "-Woverloaded-virtual"
            "-Wformat=2")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  add_options(ALL_LANGUAGES
              ALL_BUILDS
              "-Wcast-align"
              "-Wdouble-promotion"
              "-Wold-style-cast"
              "-Wnull-dereference"
              "-Wpedantic")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  add_options(ALL_LANGUAGES ALL_BUILDS -Wloop-analysis)
endif()

# FIXME TCL: compile with 2 flags below "-Wconversion" "-Wsign-conversion"

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  add_options(ALL_LANGUAGES
              ALL_BUILDS
              "-Wmisleading-indentation"
              "-Wduplicated-cond"
              "-Wlogical-op"
              "-Wuseless-cast")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
    add_options(ALL_LANGUAGES ALL_BUILDS "-Wno-implicit-fallthrough")
  endif()

endif()

set(CMAKE_CXX_FLAGS "${CXX_DIALECT_OPT_CXX11} ${CMAKE_CXX_FLAGS}")

if(TARGET_NATIVE_ARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_OPT_ARCH_NATIVE}")
endif()

# Code coverage
set(CODECOV_GCOVR_OPTIONS "-v --filter='.*/src/steps/.*'")
include(UseCodeCoverage)
if(ENABLE_CODECOVERAGE)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/${CODECOV_HTMLOUTPUTDIR}/
          DESTINATION share/steps/${STEPS_VERSION}/coverage)
  install(FILES ${CMAKE_BINARY_DIR}/${CODECOV_XMLOUTPUTFILE}
          DESTINATION share/steps/${STEPS_VERSION}/coverage)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Requirements
# -------------------------------------------------------------------
find_package(PythonInterp 3.6 REQUIRED)
find_package(PythonInterpLibs)

find_package(Threads REQUIRED)

# if USE_MPI is
#
# * False:   we do not search it
# * True:    we search and require it
# * Default: we search it, and we use it only if found
if(USE_MPI MATCHES "^(Default|True)$")
  if(USE_MPI STREQUAL "True")
    find_package(MPI REQUIRED)
  elseif(USE_MPI STREQUAL "Default")
    find_package(MPI)
  endif()
  if(MPI_FOUND)
    add_definitions(-DUSE_MPI)
    find_program(MPIRUN
                 NAMES mpirun
                 HINTS ${MPI_HOME} $ENV{MPI_HOME}
                 PATHS ${_MPI_PREFIX_PATH}
                 PATH_SUFFIXES bin
                 DOC "Executable for running MPI programs.")
  endif()
endif()

# see USE_MPI options above
if(USE_PETSC MATCHES "^(Default|True)$")
  include(FindPkgConfig)
  if(USE_PETSC STREQUAL "True")
    if(NOT MPI_FOUND)
      message(FATAL_ERROR "PETSC reqires MPI!")
    endif()
    pkg_search_module(PETSC REQUIRED PETSc)
  elseif(USE_PETSC STREQUAL "Default")
    if(NOT MPI_FOUND)
      message(
        WARNING
          "Building without PETSC. If desired enable finding MPI libs as well!")
    else()
      pkg_search_module(PETSC QUIET PETSc)
    endif()
  endif()
  if(PETSC_FOUND)
    message(STATUS "PETSc found: ${PETSC_LIBRARY_DIRS}")
    add_definitions(-DUSE_PETSC)
  endif()
endif()

# OpenBlas only made it to CMake in 3.6. Until then we include version from 3.6
find_package(BLAS REQUIRED)

# ignore LAPACK if OpenBlas is found
if(NOT BLA_VENDOR STREQUAL "OpenBLAS")
  find_package(LAPACK)
else() # if OpenBlas is used, then we have Lapack
  set(LAPACK_FOUND true)
endif()

# OpenMP
find_package(OpenMP)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  include(AppleClangPatch)
endif()
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_FLAGS "${CMAKE_FLAGS} ${OpenMP_FLAGS}")
  include_directories(SYSTEM ${OpenMP_CXX_INCLUDE_DIR})
  list(APPEND libsteps_link_libraries ${OpenMP_libomp_LIBRARY})
endif()

# Linking Libs - required in src and pysteps
list(APPEND libsteps_link_libraries ${BLAS_LIBRARIES})
if(MPI_FOUND)
  list(APPEND libsteps_link_libraries ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
  set(MPI_FOUND_HEADERS ${MPI_C_INCLUDE_PATH})
endif()
if(LAPACK_FOUND)
  list(APPEND libsteps_link_libraries ${LAPACK_LIBRARIES})
endif()
if(PETSC_FOUND)
  include_directories(SYSTEM ${PETSC_INCLUDE_DIRS})
  link_directories(${PETSC_LIBRARY_DIRS})
  link_directories(${PETSC_STATIC_LIBRARY_DIRS})
  list(APPEND libsteps_link_libraries ${PETSC_LIBRARIES})
endif()

if(NOT USE_BUNDLE_SUNDIALS)
  list(APPEND libsteps_link_libraries ${SUNDIALS_LIBRARIES})
endif(NOT USE_BUNDLE_SUNDIALS)

if(USE_BDSYSTEM_LAPACK AND NOT LAPACK_FOUND)
  message(STATUS "Unable to find LAPACK; will not build BDSystem/Lapack code.")
endif()

# Python Install Path ####
set(_default_prefixes "/usr/local" "/usr")
list(FIND _default_prefixes ${CMAKE_INSTALL_PREFIX} _index)
if(${_index} GREATER -1)
  # Default installations install Python to default site-packages
  set(_py_install_prefix ${PYTHON_SITE_PACKAGES})
else()
  set(_py_install_prefix ${CMAKE_INSTALL_PREFIX})
endif()

set(PYTHON_INSTALL_PREFIX ${_py_install_prefix}
    CACHE PATH "Path to install steps Python package")

# Python dependencies
include(FindPythonModule)
find_python_module(numpy REQUIRED)

# Python bindings tech Makes steps Python bindings based on Cython
if(PYTHONINTERP_FOUND)
  find_package(Cython 0.25.1 REQUIRED)
  add_subdirectory(pysteps)
endif()

# Makes libsteps-obj, libsteps.so
add_subdirectory(src)

# Make testing targets
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

message(
  STATUS
    "STEPS PYTHON packages will be installed to ${PYTHON_INSTALL_PREFIX}\nTo change use -DPYTHON_INSTALL_PREFIX=<dir>"
  )

# Check whether a steps version < 3.6 is already installed
find_python_module(steps)
if(STEPS_FOUND AND (STEPS_LOCATION STREQUAL "${PYTHON_INSTALL_PREFIX}/steps") AND (STEPS_VERSION_STRING VERSION_LESS 3.6.0))
    message(WARNING "STEPS ${STEPS_VERSION_STRING} is already installed in:\n${STEPS_LOCATION}\nSince STEPS 3.6.0, the layout of files in the python package has changed. It is advised to remove the previous installation before installing STEPS ${STEPS_VERSION} there. This can be done with:\n[sudo] rm -rf ${STEPS_LOCATION}")
endif()

