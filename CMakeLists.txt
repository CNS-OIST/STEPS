cmake_minimum_required(VERSION 2.8.9)
if(POLICY CMP0042)
    # Use rpath on Mac OS X
    cmake_policy(SET CMP0042 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

project(STEPS)

include(CMake/Compiler.cmake)

# C++ support in g++ pre-version 4.8 is insufficient
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    message(FATAL_ERROR "C++ compiler is g++ ${CMAKE_CXX_COMPILER_VERSION}; require version 4.8 or higher.")
endif()


# Nuke any -DNDEBUG in the compiler options introduced by CMake.
include(CMake/ManipulateVariables.cmake)
match_variables(flag_vars "^CMAKE_(C|CXX)_FLAGS.*")
foreach(var_name ${flag_vars})
    remove_word(${var_name} "-DNDEBUG")
endforeach()


# for debug symbol
include(CMake/UseDebugSymbols.cmake)

set(CMAKE_CXX_FLAGS "${CXX_DIALECT_OPT_CXX11} ${COMPILER_OPT_ARCH_NATIVE} ${CMAKE_CXX_FLAGS}")

set(VERSION_MAJOR 3)
set(VERSION_MINOR 0)
set(VERSION_PATCH 1)

set(STEPS_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# (Use modified CMake FindSWIG module to prefentially use SWIG3.0)
find_package(SWIG3 REQUIRED)
find_package(PythonInterpLibs REQUIRED)
find_package(MPI)
find_package(BLAS REQUIRED)
find_package(LAPACK)
find_package(NumPy)

if(NOT NUMPY_FOUND)
    message(STATUS "Unable to find numpy; will build STEPS without numpy support.")
endif()

if(NOT MPI_FOUND)
    message(STATUS "Unable to find MPI; will build STEPS without MPI modules.")
endif()

option(USE_BDSYSTEM_LAPACK "Use new BDSystem/Lapack code for E-Field solver" OFF)

if(USE_BDSYSTEM_LAPACK AND NOT LAPACK_FOUND)
    message(STATUS "Unable to find LAPACK; will not build BDSystem/Lapack code.")
endif()

# Makes libsteps-obj, libsteps.so
add_subdirectory(src)

# Makes _steps_swig_numpy.so or _steps_swig.so
add_subdirectory(swig)

# Makes steps python package including swig-generated code
add_subdirectory(py)

# Make testing targets
#add_subdirectory(test)
