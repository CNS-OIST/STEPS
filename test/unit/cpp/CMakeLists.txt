foreach(test_name
        # point3d
        bbox
        # tetmesh
        membership
        checkid
        # rng
        sample
        small_binomial)
  add_executable("test_${test_name}" "test_${test_name}.cpp")
  list(APPEND tests ${test_name})
endforeach()

if(BUILD_STOCHASTIC_TESTS)
    foreach(test_name
            small_binomial_stochastic)
        add_executable("test_${test_name}" "test_${test_name}.cpp")
        list(APPEND tests ${test_name})
    endforeach()
endif()

# if Lapack is used, add test for it
if(LAPACK_FOUND)
  list(APPEND tests bdsystem)
  add_executable(test_bdsystem test_bdsystem.cpp lapack_common.cpp)
endif()

if(PETSC_FOUND)
  foreach(test_name petscsystem)
    add_executable("test_${test_name}" "test_${test_name}.cpp")
    list(APPEND tests ${test_name})
  endforeach()
endif()

add_custom_target(tests)

foreach(test_name IN LISTS tests)
  string(REPLACE " "
                 ";"
                 valgrind
                 "${VALGRIND}")
  set(test_target "test_${test_name}")
  set(test_str ${valgrind} ${MPIRUN} -n 2 ${CMAKE_CURRENT_BINARY_DIR}/${test_target})

  target_link_libraries("${test_target}"
                        ${CMAKE_THREAD_LIBS_INIT}
                        ${CMAKE_DL_LIBS}
                        ${libs})

  add_test(NAME "${test_name}" COMMAND ${test_str})
  add_dependencies(tests "${test_target}")
endforeach()

# math unit test
add_executable(stepsmathTest test_point3d.cpp)
target_link_libraries(stepsmathTest stepsmath gtest_main)
add_test(NAME mathTest COMMAND stepsmathTest)

# geo unit test
add_executable(stepsgeoTest test_tetmesh.cpp)
target_link_libraries(stepsgeoTest stepsgeo gtest_main)
add_test(NAME geoTest COMMAND stepsgeoTest)

# rng unit test
add_executable(stepsrngTest test_tetmesh.cpp)
target_link_libraries(stepsrngTest stepsgeo stepsrng gtest_main)
add_test(NAME rngTest COMMAND stepsrngTest)
