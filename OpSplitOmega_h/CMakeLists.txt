add_definitions(-DPETSC_HAVE_BROKEN_RECURSIVE_MACRO)
add_executable(
  OpSplitOmega_h
  kproc/diffusions.cpp
  kproc/diffusions.hpp
  kproc/kproc_state.cpp
  kproc/kproc_state.hpp
  kproc/reactions.cpp
  kproc/reactions.hpp
  kproc/surface_reactions.hpp
  kproc/surface_reactions.cpp
  main.cpp
  mesh.cpp
  mesh.hpp
  mesh_utils.cpp
  mesh_utils.hpp
  operator/diffusion_operator.cpp
  operator/diffusion_operator.hpp
  operator/rssa_operator.cpp
  operator/rssa_operator.hpp
  operator/ssa_operator.cpp
  operator/ssa_operator.hpp
  simulation.cpp
  simulation.hpp
  simulation_data.hpp
  simulation_data.hpp
  test_rssa.cpp
  test_rssa.hpp)
target_link_libraries(OpSplitOmega_h PUBLIC zee_common Omega_h::omega_h)
install(TARGETS OpSplitOmega_h RUNTIME DESTINATION bin)

if(Zee_BUILD_TESTING)
  add_custom_command(
    OUTPUT square.msh
    COMMAND
      ${Gmsh_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/square.geo -format msh2 -2 -o square.msh
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/square.geo)
  add_custom_command(
    OUTPUT box.msh
    COMMAND
      ${Gmsh_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/box.geo -format msh2 -3 -o box.msh
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/box.geo)
  add_custom_command(
    OUTPUT square_242elts.msh
    COMMAND
      ${Gmsh_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/square.geo -setnumber elt_size 0.1 -format msh2
      -2 -o square_242elts.msh
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/square.geo)
  add_custom_command(
    OUTPUT box_749elts.msh
    COMMAND
      ${Gmsh_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/box.geo -setnumber elt_size 0.2 -format msh2 -3
      -o box_749elts.msh
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/box.geo)
  add_custom_command(
    OUTPUT cube.msh
    COMMAND
      ${Gmsh_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cube.geo -format msh2 -3 -o cube.msh
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/cube.geo)

  add_custom_target(run_OpSplitOmega_h ALL DEPENDS square.msh box.msh square_242elts.msh
                                                   box_749elts.msh cube.msh)
  add_test(NAME OpSplitSurfReacTest COMMAND $<TARGET_FILE:OpSplitOmega_h> --test 4 --scale 4.5e-6
                                            ../../OpSplitDMPlex/tests/meshes/1_tet.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_2D
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 square.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_3D
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 box.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_2D_RSSA
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 --use-rssa square.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_3D_RSSA
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 --use-rssa box.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_3D_SurfReac
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 --test -3 cube.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_2D_RSSA_CMP
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 --test -2 square_242elts.msh)
  add_mpi_test(
    NAME OpSplitOmega_h_3D_RSSA_CMP
    NUM_PROCS 2
    COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 --test -2 box_749elts.msh)
  add_test(NAME OpSplitOmega_h_3D_SREAC_2TETS
           COMMAND $<TARGET_FILE:OpSplitOmega_h> --num-iterations 100 --test 5
                   ../../OpSplitDMPlex/tests/meshes/2_tets.msh)
endif()
