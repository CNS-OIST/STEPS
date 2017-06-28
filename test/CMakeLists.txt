set(gtest_src "${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest/googletest")
set(gtest_bin "${CMAKE_CURRENT_BINARY_DIR}/third_party/gtest")

if(NOT EXISTS ${gtest_src}/CMakeLists.txt)
	message(FATAL_ERROR " Submodules do not seems to be initialized \n Please run 'git submodule update --init --recursive'")
endif()

add_subdirectory(${gtest_src} ${gtest_bin} EXCLUDE_FROM_ALL)

set(inc_dirs "${CMAKE_SOURCE_DIR}/src" "${gtest_src}/include")
set(libs libsteps gtest_main)

if (LAPACK_FOUND)
    list(APPEND libs ${LAPACK_LIBRARIES})
endif()

if (MPI_FOUND)
    list(APPEND libs ${MPI_CXX_LIBRARIES})
    list(APPEND inc_dirs ${MPI_CXX_INCLUDE_PATH})
    if (PETSC_FOUND)
        list(APPEND libs ${PETSC_LIBRARIES})
        list(APPEND inc_dirs ${PETSC_INCLUDE_DIRS})
    endif()
endif()

include_directories(${inc_dirs})

add_subdirectory(unit)
add_subdirectory(validation)
