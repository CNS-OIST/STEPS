add_library(
  stepsdist STATIC
  definition/compdef.cpp
  definition/diffdef.cpp
  definition/efield.cpp
  definition/patchdef.cpp
  definition/reacdef.cpp
  definition/sreacdef.cpp
  definition/statedef.cpp
  kproc/diffusions.cpp
  kproc/diffusions.hpp
  kproc/event_queue.cpp
  kproc/event_queue.hpp
  kproc/kproc_id.hpp
  kproc/kproc_state.cpp
  kproc/kproc_state.hpp
  kproc/propensities.cpp
  kproc/propensities.hpp
  kproc/reactions.cpp
  kproc/reactions.hpp
  kproc/reactions_iterator.hpp
  kproc/surface_reactions.cpp
  kproc/surface_reactions.hpp
  mol_state.cpp
  operator/diffusion_operator.cpp
  operator/rssa_operator.cpp
  operator/ssa_operator.cpp
  simulation.cpp
  tetopsplit.cpp)

if(PETSC_FOUND AND USE_PETSC)
  target_sources(stepsdist PRIVATE operator/efield_operator.cpp operator/efield_operator.hpp)
  target_link_libraries(stepsdist PUBLIC PETSC::PETSC)
endif()

set_property(TARGET stepsdist PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsdist PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsdist PUBLIC stepsutil stepssolver stepsgeomdist MPI::MPI_CXX
                                       Omega_h::omega_h)
