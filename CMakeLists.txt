cmake_minimum_required(VERSION 3.10.2)
if(POLICY CMP0042)
  # Use rpath on Mac OS X
  cmake_policy(SET CMP0042 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

project(STEPS LANGUAGES CXX)

set(VERSION_MAJOR 4)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(STEPS_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# hpc coding conventions
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/CMake/hpc-coding-conventions/cpp")
  message(FATAL_ERROR "Submodules do not seem to be initialized. Please run 'git submodule update --init --recursive'")
endif()
set(${PROJECT_NAME}_ClangFormat_EXCLUDES_RE
    ".*/pysteps/.*$$"
    ".*/src/third_party/.*$$"
    ".*/test/third_party/.*$$"
    CACHE STRING "list of regular expressions to exclude C/C++ files from formatting" FORCE)
set(${PROJECT_NAME}_ClangTidy_EXCLUDES_RE
    ".*/pysteps/.*$$"
    ".*/src/third_party/.*$$"
    ".*/test/third_party/.*$$"
    CACHE STRING "list of regular expressions to exclude C/C++ files from formatting" FORCE)
add_subdirectory(CMake/hpc-coding-conventions/cpp)

# OPTIONS
option(BUILD_STOCHASTIC_TESTS "Build stochastic tests" ON)

set(TARGET_NATIVE_ARCH_DEFAULT ON)
if(APPLE AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64") AND (CMAKE_CXX_COMPILER_ID STREQUAL AppleClang))
  # Disable native instructions on Apple M1 with AppleClang since it is not supported.
  set(TARGET_NATIVE_ARCH_DEFAULT OFF)
endif()
option(TARGET_NATIVE_ARCH "Generate non-portable arch-specific code" ${TARGET_NATIVE_ARCH_DEFAULT})
option(USE_BDSYSTEM_LAPACK "Use new BDSystem/Lapack code for E-Field solver"
       OFF)
option(STEPS_USE_DIST_MESH "Add solvers based on distributed mesh" ON)
option(USE_64_BITS_INDICES "Use 64bits indices instead of 32" OFF)
option(STEPS_ENABLE_ERROR_ON_WARNING "Add -Werror to STEPS compilation" OFF)

option(USE_BUNDLE_EASYLOGGINGPP "Use bundled version of easylogging" ON)
option(USE_BUNDLE_RANDOM123 "Use bundled version of random123" ON)
option(USE_BUNDLE_SUNDIALS "Use bundled version of cvode" ON)
option(USE_BUNDLE_OMEGA_H "Use bundled version of Omega_h" ON)

option(STEPS_USE_CALIPER_PROFILING "Use Caliper instrumentation" OFF)
option(STEPS_USE_LIKWID_PROFILING "Use Likwid instrumentation" OFF)
option(STEPS_USE_NATIVE_PROFILING "Use STEPS region tracker" OFF)
# mark as advanced profiling options
include(CMake/MarkAsAdvancedAll.cmake)
mark_as_advanced_all("PROFILING")

set(USE_MPI "Default" CACHE STRING "Use MPI for parallel solvers")
# Toggleable options
set_property(CACHE USE_MPI PROPERTY STRINGS "Default;True;False")
set(USE_PETSC "Default"
    CACHE STRING "Use PETSC library for parallel E-Field solver")
set_property(CACHE USE_PETSC PROPERTY STRINGS "Default;True;False")

if(USE_BUNDLE_EASYLOGGINGPP)
  set(EASYLOGGINGPP_INCLUDE_DIRS
      "${CMAKE_SOURCE_DIR}/src/third_party/easyloggingpp/src")
  # mark as advanced gtest stuff
  include(CMake/MarkAsAdvancedAll.cmake)
  mark_as_advanced_all("gtest")
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

if(STEPS_USE_DIST_MESH)
  if(USE_BUNDLE_OMEGA_H)
    include(BundleOmega_h)
  else()
    find_package(Omega_h REQUIRED)
    string(REPLACE "${Omega_h_CXX_FLAGS}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()
  add_definitions(-DSTEPS_USE_DIST_MESH)
endif()

if(USE_BUNDLE_SUNDIALS)
  include(BundleSUNDIALS)
else()
  set(LOAD_CVODE ON)
  find_package(SUNDIALS 2.6...3.9.9 REQUIRED)
  if(CMAKE_VERSION VERSION_LESS 3.19.0)
    # For older versions of cmake we need to check upper limit manually
    if(SUNDIALS_VERSION VERSION_GREATER_EQUAL 4.0.0)
      message(SEND_ERROR "SUNDIALS is supported up to version 3.x but found ${SUNDIALS_VERSION}")
    endif()
  endif()
  include(ImportModernLib)
  add_library_target(NAME SUNDIALS
    LIBRARIES ${SUNDIALS_LIBRARIES}
    INCLUDE_DIRECTORIES ${SUNDIALS_INCLUDE_DIR})
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
  add_options(ALL_LANGUAGES ALL_BUILDS "-Wno-deprecated-declarations")
  add_options(ALL_LANGUAGES ALL_BUILDS "-Wno-format-nonliteral")
  add_options(ALL_LANGUAGES ALL_BUILDS -Wloop-analysis)
  add_options(ALL_LANGUAGES ALL_BUILDS "-Wno-unused-command-line-argument")
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

set(CMAKE_CXX_FLAGS "${CXX_DIALECT_OPT_CXX17} ${CMAKE_CXX_FLAGS}")

if(TARGET_NATIVE_ARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_OPT_ARCH_NATIVE}")
endif()

include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/src/third_party/gsl-lite/include)

# Code coverage
set(CODECOV_GCOVR_OPTIONS "-v --filter='.*/src/steps/.*'")
include(UseCodeCoverage)
if(ENABLE_CODECOVERAGE)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/${CODECOV_HTMLOUTPUTDIR}/
          DESTINATION share/steps/${STEPS_VERSION}/coverage)
  install(FILES ${CMAKE_BINARY_DIR}/${CODECOV_XMLOUTPUTFILE}
          DESTINATION share/steps/${STEPS_VERSION}/coverage)
endif()

set(CMAKE_C_FLAGS_RELWITHDEBINFOANDASSERT "-O2 -g" CACHE STRING
  "C compiler flags corresponding to RelWithDebInfo build type with assert enabled." FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFOANDASSERT "-O2 -g" CACHE STRING
  "C++ compiler flags corresponding to RelWithDebInfo build type with assert enabled." FORCE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Toggleable build type
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel;RelWithDebInfoAndAssert")

# Requirements
# -------------------------------------------------------------------
find_package(PythonInterp 3.6 REQUIRED)
find_package(PythonInterpLibs)


if(APPLE)
  # assume built-in pthreads on MacOS
  # https://stackoverflow.com/questions/54587052
  set(CMAKE_THREAD_LIBS_INIT "-lpthread")
  set(CMAKE_HAVE_THREADS_LIBRARY 1)
  set(CMAKE_USE_WIN32_THREADS_INIT 0)
  set(CMAKE_USE_PTHREADS_INIT 1)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
else()
    find_package(Threads REQUIRED)
endif()

find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

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
    # this hint tries to ensure that the mpirun selected is compatible with the compiler
    get_filename_component(MPI_BIN_PATH ${MPI_CXX_COMPILER} DIRECTORY)
    find_program(MPIRUN
                 NAMES srun mpirun
                 HINTS ${MPI_BIN_PATH} ${MPI_HOME} $ENV{MPI_HOME}
                 PATHS ${_MPI_PREFIX_PATH}
                 PATH_SUFFIXES bin
                 DOC "Executable for running MPI programs.")
    if(MPI_CXX_COMPILE_OPTIONS STREQUAL "-framework")
      # reset incomplete option added on MacOS in some cases
      set(MPI_CXX_COMPILE_OPTIONS "" CACHE STRING "" FORCE)
    endif()
  endif()
endif()

# see USE_MPI options above
if(USE_PETSC MATCHES "^(Default|True)$")
  include(FindPkgConfig)
  if(USE_PETSC STREQUAL "True")
    if(NOT MPI_FOUND)
      message(FATAL_ERROR "PETSC requires MPI!")
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

if(STEPS_USE_DIST_MESH AND USE_PETSC)
  if(PETSC_SIZEOF_PETSCINT EQUAL 4)
    message(WARNING
    "Distributed mesh requires PETSc to be compiled with configure option "
    "'--with-64-bit-indices=1' to support mesh with "
    "more than 2.1 billions tetrahedrons. This option is recommended anyway "
    "on clusters for efficiency purpose but not required for smaller meshes.")
  endif()
endif()

# OpenBlas only made it to CMake in 3.6. Until then we include version from 3.6
find_package(BLAS REQUIRED)
# Look for lapack
find_package(LAPACK)

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

# Libs - required in src and pysteps
include(ImportModernLib)
add_library_target(NAME BLAS
  LIBRARIES ${BLAS_LIBRARIES})

if(LAPACK_FOUND)
  add_library_target(NAME LAPACK
    LIBRARIES ${LAPACK_LIBRARIES})
endif()

if(PETSC_FOUND)
  # needed before cmake 3.12
  if(NOT PETSC_LINK_LIBRARIES)
    set(PETSC_LINK_LIBRARIES "${PETSC_LIBRARY_DIRS}/lib${PETSC_LIBRARIES}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  endif()
  add_library_target(NAME PETSC
    LIBRARIES ${PETSC_LINK_LIBRARIES}
    INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
endif()

if(USE_BDSYSTEM_LAPACK AND NOT LAPACK_FOUND)
  message(STATUS "Unable to find LAPACK; will not build BDSystem/Lapack code.")
endif()

# Python Install Path ####
set(_default_prefixes "/usr/local" "/usr")
list(FIND _default_prefixes ${CMAKE_INSTALL_PREFIX} _index)
if(${_index} GREATER -1)
  # Default installations install Python to default site-packages
  set(_py_install_prefix ${PYTHON_SITE_PACKAGES})
  get_filename_component(bundle_lib_install_dir "${PYTHON_LIBRARY}" DIRECTORY)
else()
  set(_py_install_prefix ${CMAKE_INSTALL_PREFIX})
  set(bundle_lib_install_dir ${CMAKE_INSTALL_PREFIX}/lib)
endif()

set(PYTHON_INSTALL_PREFIX ${_py_install_prefix}
    CACHE PATH "Path to install steps Python package")

# Python dependencies
include(FindPythonModule)
find_python_module(numpy REQUIRED)
if(MPI_FOUND AND STEPS_USE_DIST_MESH)
  find_python_module(mpi4py REQUIRED)
endif()

# Python bindings tech Makes steps Python bindings based on Cython
if(PYTHONINTERP_FOUND)
  find_package(Cython 0.25.1 REQUIRED)
  add_subdirectory(pysteps)
endif()

# Install the bundle libraries if any
if(USE_BUNDLE_OMEGA_H)
  # Add transitive link of Gmsh when linking against Omega_h
  find_package(Gmsh REQUIRED)
  include_directories("${Gmsh_INCLUDE_DIRS}")
  install(FILES "${tpdir}/omega_h-install/lib/libomega_h${CMAKE_SHARED_LIBRARY_SUFFIX}" DESTINATION "${bundle_lib_install_dir}")
endif()
if(USE_BUNDLE_SUNDIALS)
  install(DIRECTORY "${tpdir}/SUNDIALS-install/lib/" DESTINATION "${bundle_lib_install_dir}")
endif()
unset(python_library_dir)


# Makes libsteps-obj, libsteps.so
add_subdirectory(src)

# Make testing targets
include(CTest)
include(AddPythonTests)
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
