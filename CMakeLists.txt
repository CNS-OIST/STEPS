cmake_minimum_required(VERSION 3.6)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

project(
  Zee
  VERSION 0.1.0
  LANGUAGES CXX)

include(FindPkgConfig)
option(Zee_BUILD_TESTING "Build tests" ON)

if(Zee_BUILD_TESTING)
  enable_testing()
endif()

add_subdirectory(hpc-coding-conventions/cpp)

bob_begin_package()

bob_option(Zee_USE_PETSc "Build examples using PETSc" ON)

bob_option(Zee_USE_MPI "Enable MPI" ON)
if(Zee_USE_MPI)
  find_package(MPI REQUIRED)
  find_program(
    MPIRUN
    NAMES srun mpirun
    HINTS ${MPI_HOME} $ENV{MPI_HOME}
    PATHS ${_MPI_PREFIX_PATH}
    PATH_SUFFIXES bin
    DOC "Executable for running MPI programs.")

  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif()

find_package(OpenMP)

set(Zee_USE_Omega_h_DEFAULT ON)
bob_public_dep(Omega_h)

find_package(Boost REQUIRED COMPONENTS timer)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

pkg_search_module(PETSC REQUIRED PETSc)
include_directories(SYSTEM ${PETSC_INCLUDE_DIRS})
link_directories(${PETSC_LIBRARY_DIRS})
link_directories(${PETSC_STATIC_LIBRARY_DIRS})

bob_option(Zee_USE_TIMEMORY "Enable timemory API" OFF)
if(Zee_USE_TIMEMORY)
  find_package(timemory REQUIRED)
  include_directories(SYSTEM ${timemory_INCLUDE_DIR})
  link_directories(${timemory_LIB_DIR})
  add_definitions(-DZee_USE_TIMEMORY)
endif()

include_directories(SYSTEM third_party)
include_directories(SYSTEM third_party/hadoken/include)
include_directories(common/include)
include_directories(${CMAKE_BINARY_DIR}/common/include)

bob_begin_cxx_flags()
bob_cxx14_flags()

if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_FLAGS "${CMAKE_FLAGS} ${OpenMP_FLAGS}")
else()
  message(WARNING "OpenMP support is disabled because it could not be found.")
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-weak-template-vtables")
endif()
bob_end_cxx_flags()

function(add_mpi_test)
  set(options "")
  set(oneValueArgs NAME NUM_PROCS)
  set(multiValueArgs COMMAND CONFIGURATIONS)
  cmake_parse_arguments(add_mpi_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT add_mpi_test_NUM_PROCS)
    set(add_mpi_test_NUM_PROCS 1)
  endif()
  if(MPI_FOUND)
    list(GET add_mpi_test_COMMAND 0 add_mpi_test_executable)
    list(REMOVE_AT add_mpi_test_COMMAND 0)
    set(TEST_STR
        ${MPIRUN}
        ${MPIEXEC_NUMPROC_FLAG}
        ${add_mpi_test_NUM_PROCS}
        ${MPIEXEC_PREFLAGS}
        ${add_mpi_test_executable}
        ${MPIEXEC_POSTFLAGS}
        ${add_mpi_test_COMMAND})
  else(MPI_FOUND)
    if(NOT (${add_mpi_test_NUM_PROCS} EQUAL "1"))
      message(STATUS "test ${add_mpi_test_NAME} ignored because MPIRUN not found!")
      return()
    endif()
    set(TEST_STR ${add_mpi_test_COMMAND})
  endif(MPI_FOUND)
  set(add_test_args NAME ${add_mpi_test_NAME} COMMAND ${TEST_STR})
  if(add_mpi_test_CONFIGURATIONS)
    list(APPEND add_test_args CONFIGURATIONS ${add_mpi_test_CONFIGURATIONS})
  endif()
  if(add_mpi_test_WORKING_DIRECTORY)
    list(APPEND add_test_args WORKING_DIRECTORY ${add_mpi_test_WORKING_DIRECTORY})
  endif()
  add_test(${add_test_args})
endfunction(add_mpi_test)

bob_get_commit()

add_subdirectory(common)
add_subdirectory(OpSplitOmega_h)
add_subdirectory(test)
bob_end_package()
