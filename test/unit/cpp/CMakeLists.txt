if(APPLE)
  # so that built executables here can find the bundled shared libraries copied there.
  list(APPEND CMAKE_BUILD_RPATH ${CMAKE_BINARY_DIR}/lib/steps)
endif(APPLE)

function(test_unit)
  set(multiValueArgs TARGETS DEPENDENCIES INCLUDE MPI_RANKS)
  set(oneValueArgs PATH)
  cmake_parse_arguments(TEST_UNIT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT TEST_UNIT_PATH)
    set(TEST_UNIT_PATH ".")
  endif()

  string(REPLACE " " ";" valgrind "${VALGRIND}")
  foreach(target ${TEST_UNIT_TARGETS})
    if(TEST_UNIT_MPI_RANKS)
      foreach(ranks ${TEST_UNIT_MPI_RANKS})
        set(test_str ${valgrind} ${MPIRUN} -n ${ranks}
                     "${CMAKE_CURRENT_BINARY_DIR}/test_${target}_${ranks}")
        add_executable("test_${target}_${ranks}" "${TEST_UNIT_PATH}/test_${target}.cpp")
        if(TEST_UNIT_DEPENDENCIES)
          target_link_libraries("test_${target}_${ranks}" "${TEST_UNIT_DEPENDENCIES}")
        endif()
        if(TEST_UNIT_INCLUDE)
          target_include_directories("test_${target}_${ranks}" PRIVATE "${TEST_UNIT_INCLUDE}")
        endif()
        add_test(NAME "${target}_${ranks}" COMMAND ${test_str})
        set_tests_properties("${target}_${ranks}" PROPERTIES PROCESSORS ${ranks})
        if(STEPS_SANITIZERS)
          cpp_cc_configure_sanitizers(TEST "${target}_${ranks}" PRELOAD)
        endif(STEPS_SANITIZERS)
      endforeach()
    else()
      set(test_str ${valgrind} "${CMAKE_CURRENT_BINARY_DIR}/test_${target}")
      add_executable("test_${target}" "${TEST_UNIT_PATH}/test_${target}.cpp")
      if(TEST_UNIT_DEPENDENCIES)
        target_link_libraries("test_${target}" "${TEST_UNIT_DEPENDENCIES}")
      endif()
      if(TEST_UNIT_INCLUDE)
        target_include_directories("test_${target}" PRIVATE "${TEST_UNIT_INCLUDE}")
      endif()
      add_test(NAME "${target}" COMMAND ${test_str})
      if(STEPS_SANITIZERS)
        cpp_cc_configure_sanitizers(TEST "${target}" PRELOAD)
        message(STATUS "Configure sanitizers for test: ${target}")
      endif(STEPS_SANITIZERS)
    endif()
  endforeach()
endfunction()

test_unit(
  TARGETS checkid
          checkpoint
          errorlog
          flat_multimap
          membership
          memory_tracker
          optional_num
          strongid
          time_tracker
  PATH util
  DEPENDENCIES stepsutil gtest_main)

test_unit(TARGETS bbox point3d sample ghkcurr DEPENDENCIES stepsmath gtest_main)

test_unit(TARGETS tetmesh DEPENDENCIES stepsgeom gtest_main ${METIS_LIBRARIES})

test_unit(TARGETS rng small_binomial DEPENDENCIES stepsrng gtest_main)

test_unit(TARGETS efield_tetmesh DEPENDENCIES stepssolver gtest_main)

if(BUILD_STOCHASTIC_TESTS)
  test_unit(TARGETS small_binomial_stochastic DEPENDENCIES stepsrng gtest_main)
endif()

if(PETSC_FOUND AND MPI_FOUND)
  test_unit(
    TARGETS petscsystem
    DEPENDENCIES PETSC::PETSC MPI::MPI_CXX gtest_main
    MPI_RANKS 2 3)
endif()

if(LAPACK_FOUND)
  add_library(lapack_common STATIC lapack_common.cpp)
  test_unit(TARGETS bdsystem DEPENDENCIES lapack_common stepsutil stepssolver gtest_main)
endif()

if(STEPS_USE_DIST_MESH)
  configure_file(test_common.hpp.in test_common.hpp @ONLY)

  test_unit(TARGETS propensities DEPENDENCIES stepsdist)

  test_unit(TARGETS mol_state DEPENDENCIES stepsdist gtest_main)

  test_unit(TARGETS kproc_event_queue DEPENDENCIES stepsdist gtest_main)

  include_directories(${CMAKE_CURRENT_BINARY_DIR})
endif()
