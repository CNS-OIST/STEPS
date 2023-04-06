function(test_unit)
    set(multiValueArgs TARGETS DEPENDENCIES INCLUDE MPI_RANKS)
    set(oneValueArgs PATH)
    cmake_parse_arguments(TEST_UNIT "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})

    if(NOT TEST_UNIT_PATH)
      set(TEST_UNIT_PATH ".")
    endif()

    string(REPLACE " "
                 ";"
                 valgrind
                 "${VALGRIND}")
    foreach(target ${TEST_UNIT_TARGETS})
        if(TEST_UNIT_MPI_RANKS)
            foreach(ranks ${TEST_UNIT_MPI_RANKS})
              set(test_str ${valgrind} ${MPIRUN} -n ${ranks} "${CMAKE_CURRENT_BINARY_DIR}/test_${target}_${ranks}")
              add_executable("test_${target}_${ranks}" "${TEST_UNIT_PATH}/test_${target}.cpp")
              if(TEST_UNIT_DEPENDENCIES)
                target_link_libraries("test_${target}_${ranks}" "${TEST_UNIT_DEPENDENCIES}")
              endif()
              if(TEST_UNIT_INCLUDE)
                target_include_directories("test_${target}_${ranks}" PRIVATE "${TEST_UNIT_INCLUDE}")
              endif()
              add_test(NAME "${target}_${ranks}" COMMAND ${test_str})
              set_tests_properties("${target}_${ranks}" PROPERTIES PROCESSORS ${ranks})
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
        endif()
    endforeach()
endfunction()

test_unit(TARGETS checkid
                  errorlog
                  flat_multimap
                  membership
                  memory_tracker
                  optional_num
                  peak_rss
                  strongid
                  time_tracker
          PATH util
          DEPENDENCIES stepsutil
                         gtest_main)

test_unit(TARGETS bbox
                  point3d
                  sample
                  ghkcurr
          DEPENDENCIES stepsmath
                         gtest_main)

test_unit(TARGETS tetmesh
          DEPENDENCIES stepsgeom
                         gtest_main)

test_unit(TARGETS rng
                  small_binomial
          DEPENDENCIES stepsrng
                         gtest_main)

test_unit(TARGETS efield_tetmesh DEPENDENCIES stepssolver gtest_main)

if(BUILD_STOCHASTIC_TESTS)
    test_unit(TARGETS small_binomial_stochastic
              DEPENDENCIES stepsrng
                             gtest_main)
endif()

if(PETSC_FOUND AND MPI_FOUND)
    test_unit(TARGETS petscsystem
              DEPENDENCIES PETSC::PETSC
                             MPI::MPI_CXX
                             gtest_main
              MPI_RANKS 2 3)
endif()

if(LAPACK_FOUND)
  add_library(lapack_common STATIC lapack_common.cpp)
  test_unit(TARGETS bdsystem
            DEPENDENCIES lapack_common stepsutil stepssolver gtest_main)
endif()

if(STEPS_USE_DIST_MESH)
  configure_file(test_common.hpp.in test_common.hpp @ONLY)
  test_unit(TARGETS opsplit
            DEPENDENCIES stepsdist
                           stepsdisttest
            INCLUDE "${CMAKE_CURRENT_BINARY_DIR}"
            MPI_RANKS 1 4)

  test_unit(TARGETS propensities
          DEPENDENCIES stepsdist stepsdisttest)

  test_unit(TARGETS mol_state
          DEPENDENCIES stepsdist gtest_main)

  test_unit(TARGETS kproc_event_queue
          DEPENDENCIES stepsdist gtest_main)

  include_directories(${CMAKE_CURRENT_BINARY_DIR})
  # The function above can be extended to include also args and cover
  # the dist steps test -> future work!
  function(add_mpi_test)
    set(options "")
    set(oneValueArgs NAME NUM_PROCS)
    set(multiValueArgs COMMAND CONFIGURATIONS)
    cmake_parse_arguments(add_mpi_test "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})
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
        message(
          STATUS "test ${add_mpi_test_NAME} ignored because MPIRUN not found!")
        return()
      endif()
      set(TEST_STR ${add_mpi_test_COMMAND})
    endif(MPI_FOUND)
    set(add_test_args NAME ${add_mpi_test_NAME} COMMAND ${TEST_STR})
    if(add_mpi_test_CONFIGURATIONS)
      list(APPEND add_test_args CONFIGURATIONS ${add_mpi_test_CONFIGURATIONS})
    endif()
    if(add_mpi_test_WORKING_DIRECTORY)
      list(APPEND add_test_args WORKING_DIRECTORY
                  ${add_mpi_test_WORKING_DIRECTORY})
    endif()
    add_test(${add_test_args})
    set_tests_properties(${add_mpi_test_NAME} PROPERTIES PROCESSORS ${add_mpi_test_NUM_PROCS})
  endfunction(add_mpi_test)

  add_mpi_test(
    NAME
    tetopsplit_box
    NUM_PROCS
    2
    COMMAND
    $<TARGET_FILE:tetopsplit_dist>
    --num-iterations
    100
    --test
    1
    ${CMAKE_SOURCE_DIR}/test/mesh/cube.msh)
   add_mpi_test(
     NAME
     tetopsplit_SReac_Validation_MultiTet
     NUM_PROCS
     4
     COMMAND
     $<TARGET_FILE:tetopsplit_dist>
     --test
     4
     --scale
     4.5e-6
     ${CMAKE_SOURCE_DIR}/test/mesh/cube.msh)
   add_mpi_test(
     NAME
     tetopsplit_SReac_Validation_MultiTet_Indep_KProcs
     NUM_PROCS
     4
     COMMAND
     $<TARGET_FILE:tetopsplit_dist>
     --independent-kprocs
     --test
     4
     --scale
     4.5e-6
     ${CMAKE_SOURCE_DIR}/test/mesh/cube.msh)
   add_mpi_test(
     NAME
     tetopsplit_SReac_Validation_MultiTet_Gibson_Bruck
     NUM_PROCS
     4
     COMMAND
     $<TARGET_FILE:tetopsplit_dist>
     --use-gibson-bruck
     --test
     4
     --scale
     4.5e-6
     ${CMAKE_SOURCE_DIR}/test/mesh/cube.msh)
   add_mpi_test(
     NAME
     tetopsplit_SReac_Validation_MultiTet_Indep_KProcs_RSSA
     NUM_PROCS
     1
     COMMAND
     $<TARGET_FILE:tetopsplit_dist>
     --independent-kprocs
     --use-rssa
     --test
     4
     --scale
     4.5e-6
     ${CMAKE_SOURCE_DIR}/test/mesh/cube.msh)
  add_mpi_test(
    NAME tetopsplit_SREAC_2TETS
    NUM_PROCS 1
    COMMAND $<TARGET_FILE:tetopsplit_dist> --num-iterations 100 --test 5 --end-time 20.0
            ${CMAKE_SOURCE_DIR}/test/mesh/2_tets.msh)
  add_mpi_test(
    NAME tetopsplit_GHKSurfaceReactions
    NUM_PROCS 1
    COMMAND $<TARGET_FILE:tetopsplit_dist> --test 8
            ${CMAKE_SOURCE_DIR}/test/mesh/2_tets.msh)
  add_mpi_test(
          NAME tetopsplit_CaBurstIntegrationTest_1
          NUM_PROCS 1
          COMMAND $<TARGET_FILE:tetopsplit_dist> --test 13  -ksp_type pipecg -pc_type jacobi -ksp_rtol 1e-8 --scale
          1e-6 --rng-seed $i --efield-dt 1e-5 --end-time 0.001
          ${CMAKE_SOURCE_DIR}/test/mesh/3tets_2patches_2comp.msh)
  add_mpi_test(
          NAME tetopsplit_CaBurstIntegrationTest_2
          NUM_PROCS 2
          COMMAND $<TARGET_FILE:tetopsplit_dist> --test 13  -ksp_type pipecg -pc_type jacobi -ksp_rtol 1e-8 --scale
          1e-6 --rng-seed $i --efield-dt 1e-5 --end-time 0.001
          ${CMAKE_SOURCE_DIR}/test/mesh/3tets_2patches_2comp_split2/3tets_2patches_2comp)
  add_mpi_test(
    NAME testopsplit_DiffusionBoundary
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:tetopsplit_dist> --test 9
            ${CMAKE_SOURCE_DIR}/test/mesh/3_tets.msh)
  add_mpi_test(
          NAME testopsplit_InitMolDistribution
          NUM_PROCS 2
          COMMAND $<TARGET_FILE:tetopsplit_dist> --test 11 --rng-seed 3
          ${CMAKE_SOURCE_DIR}/test/mesh/diamond.msh)
endif()
