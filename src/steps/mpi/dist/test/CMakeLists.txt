add_library(stepsdisttest STATIC
    ca_burst_background.cpp
    ca_burst_full.cpp
    ca_burst_integration_test.cpp
    diff_sim.cpp
    ghk_current_unittest.cpp
    multi_comp.cpp
    multi_comp_diff.cpp
    rallpack3.cpp
    scenario.cpp
    simdef.cpp
    simple_model.cpp
    simulation.cpp
    single_comp_dist.cpp
    sreac_unittest.cpp
    sreac_validation.cpp
    validation.cpp
)

set_property(TARGET stepsdisttest PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsdisttest PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsdisttest PUBLIC stepsutil MPI::MPI_CXX Omega_h::omega_h)

if(PETSC_FOUND)
    target_link_libraries(stepsdisttest PUBLIC PETSC::PETSC)
endif()

