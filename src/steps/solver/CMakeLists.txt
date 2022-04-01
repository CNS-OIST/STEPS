add_library(stepssolver STATIC
    api_comp.cpp
    api_main.cpp
    api_memb.cpp
    api_patch.cpp
    api_tet.cpp
    api_vert.cpp
    api_tri.cpp
    api_diffboundary.cpp
    api_recording.cpp
    api_batchdata.cpp
    api_roidata.cpp
    compdef.cpp
    diffdef.cpp
    patchdef.cpp
    api_sdiffboundary.cpp
    reacdef.cpp
    specdef.cpp
    sreacdef.cpp
    statedef.cpp
    chandef.cpp
    ghkcurrdef.cpp
    diffboundarydef.cpp
    ohmiccurrdef.cpp
    vdeptransdef.cpp
    vdepsreacdef.cpp
    sdiffboundarydef.cpp
    efield/dVsolver.cpp
    efield/bdsystem.cpp
    efield/dVsolver.cpp
    efield/efield.cpp
    efield/matrix.cpp
    efield/tetcoupler.cpp
    efield/tetmesh.cpp
    efield/vertexconnection.cpp
    efield/vertexelement.cpp
)

if(LAPACK_FOUND)
    target_sources(stepssolver PRIVATE
        "efield/bdsystem_lapack.cpp"
    )
    target_link_libraries(stepssolver PUBLIC LAPACK::LAPACK BLAS::BLAS)
endif()

if(PETSC_FOUND)
    target_sources(stepssolver PRIVATE
        "efield/dVsolver_petsc.cpp"
    )
    target_link_libraries(stepssolver PUBLIC MPI::MPI_CXX PETSC::PETSC)
endif()

set_property(TARGET stepssolver PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepssolver PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepssolver PUBLIC stepsutil stepsgeom stepsmodel stepsrng)