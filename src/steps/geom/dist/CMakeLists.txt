add_library(stepsgeomdist STATIC
    distcomp.cpp
    distmemb.cpp
    distmesh.cpp
    distpatch.cpp
    measure.cpp
)

set_property(TARGET stepsgeomdist PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsgeomdist PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsgeomdist PUBLIC stepsmath stepsutil stepsmpi stepsmodel MPI::MPI_CXX)

if(PETSC_FOUND)
    target_link_libraries(stepsgeomdist PUBLIC PETSC::PETSC)
endif()
