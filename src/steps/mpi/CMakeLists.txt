add_library(stepsmpi STATIC mpi_init.cpp mpi_finish.cpp)
target_link_libraries(stepsmpi PUBLIC stepsutil MPI::MPI_CXX)

if(PETSC_FOUND)
  target_link_libraries(stepsmpi PUBLIC PETSC::PETSC)
endif()

set_property(TARGET stepsmpi PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory(tetopsplit)
add_subdirectory(tetvesicle)

if(STEPS_USE_DIST_MESH)
  add_subdirectory(dist/tetopsplit)
endif()
