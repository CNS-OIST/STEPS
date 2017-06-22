set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS "-g ${CXX_DIALECT_OPT_CXX11} -O0")

foreach(test_name point3d bbox tetmesh membership checkid rng sample small_binomial)
    add_executable("test_${test_name}" "test_${test_name}.cpp")
    list(APPEND tests ${test_name})
endforeach()

# if Lapack is used, add test for it
if(LAPACK_FOUND)
    list(APPEND tests bdsystem)
    add_executable(test_bdsystem test_bdsystem.cpp lapack_common.cpp)
    if(MPI_FOUND)
        list(APPEND tests slusystem)
        add_executable(test_slusystem test_slusystem.cpp lapack_common.cpp)
    endif()
endif()

if (PETSC_FOUND)
    foreach (test_name petscsystem)
        add_executable("test_${test_name}" "test_${test_name}.cpp")
        list(APPEND tests ${test_name})
    endforeach()
endif()


add_custom_target(tests)

foreach(test_name IN LISTS tests)
    set(test_target "test_${test_name}")

    target_link_libraries("${test_target}" ${CMAKE_THREAD_LIBS_INIT}  ${libs})

    add_test(NAME "${test_name}" COMMAND "${test_target}")
    add_dependencies(tests "${test_target}")
endforeach()

