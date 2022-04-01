# ------------------------------------------------------------------
# Definitions
# ------------------------------------------------------------------
set(steps_include_dir "${PROJECT_SOURCE_DIR}/src/steps")
set(setup_py "${CMAKE_CURRENT_BINARY_DIR}/cmake_setup.py")
set(libdir "${CMAKE_BINARY_DIR}/lib")
set(cymodname cysteps)
if(MPI_FOUND)
  if(STEPS_USE_DIST_MESH)
    set(cymodname ${cymodname}_dist)
  else()
    set(cymodname ${cymodname}_mpi)
  endif()
endif()
set(cython_outputs "${CMAKE_CURRENT_SOURCE_DIR}/steps/${cymodname}.cpp")
set(steps_package_timestamp "${CMAKE_CURRENT_BINARY_DIR}/python_built.tmp")
set(cython_src ${cymodname}.pyx)
# Track cython sources
set(cy_srcs
    cysteps.pyx
    cysteps__globals.pyx
    cysteps_geom.pyx
    cysteps_model.pyx
    cysteps_mpi.pyx
    cysteps_rng.pyx
    cysteps_solver.pyx
    std.pxd
    std_fstream.pxd
    steps.pxd
    steps_common.pxd
    steps_model.pxd
    steps_mpi.pxd
    steps_rng.pxd
    steps_solver.pxd
    steps_tetexact.pxd
    steps_tetmesh.pxd
    steps_tetode.pxd
    steps_wm.pxd
    steps_wmdirect.pxd
    steps_wmrk4.pxd
    steps_wmrssa.pxd)

if(STEPS_USE_DIST_MESH)
  list(APPEND cy_srcs cysteps_dist.pyx)
endif()

# Configure setup.py to include version
configure_file(cmake_setup.py.in ${setup_py})

# ------------------------------------------------------------------
# Make rules
# ------------------------------------------------------------------
# .pyx -> .cpp
if(PETSC_FOUND AND USE_PETSC)
  set(cython_use_petsc True)
else()
  set(cython_use_petsc False)
endif()
add_custom_command(
  OUTPUT ${cython_outputs}
  COMMAND
    ${CYTHON_EXECUTABLE} --fast-fail --cplus -I ${steps_include_dir}
    --compile-time-env USE_PETSC=${cython_use_petsc}
    --compile-time-env STEPS_USE_DIST_MESH=${STEPS_USE_DIST_MESH}
    --output-file ${cython_outputs} ${cython_src} -X embedsignature=True
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS libsteps ${cy_srcs}
  COMMENT "Cythonizing extension ${cython_src}")

# .cpp -> .so (it will be a single .so without dependencies on libsteps)
set(CYSTEPS_COMPILE_FLAGS
    "-Wno-strict-aliasing -Wno-tautological-compare \
     -Wno-format -Wno-double-promotion -Wno-unsequenced \
     -Wno-format-nonliteral \
     -Wno-unused-parameter -Wno-old-style-cast -Wno-unused-variable \
     -Wno-shadow -Wno-parentheses -Wno-cast-align")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CYSTEPS_COMPILE_FLAGS "${CYSTEPS_COMPILE_FLAGS} -Wno-useless-cast -Wno-deprecated-declarations -Wno-overflow")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
    set(CYSTEPS_COMPILE_FLAGS "${CYSTEPS_COMPILE_FLAGS} -Wno-cast-function-type")
  endif()
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  set(CYSTEPS_COMPILE_FLAGS "${CYSTEPS_COMPILE_FLAGS} -Wno-c++17-extensions -Wno-implicit-fallthrough")
endif()

python_add_module(cysteps ${cython_outputs})

set_source_files_properties(${cython_outputs} PROPERTIES COMPILE_FLAGS
                            "${CYSTEPS_COMPILE_FLAGS}")

target_include_directories(cysteps
                           PRIVATE ${PYTHON_INCLUDE_DIRS}
                                   ${steps_include_dir})
target_link_libraries(cysteps libsteps_static ${libsteps_link_libraries})

set_target_properties(
  cysteps PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${libdir}/steps" OUTPUT_NAME
                     ${cymodname})

# Pure python package
add_custom_target(steps_package_pure_py ALL DEPENDS ${steps_package_timestamp}
)# random file bridging target->command
add_custom_command(
  OUTPUT ${steps_package_timestamp}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND
    ${PYTHON_EXECUTABLE} ${setup_py} build --build-base ${CMAKE_BINARY_DIR}
    install --install-lib "${libdir}" --no-compile --force
  COMMAND
    ${CMAKE_COMMAND} -E touch ${steps_package_timestamp}





    # setup.py already tracks changed .py files. We set DEPENDS just to set the
    # order: copy .py after building the lib
  DEPENDS cysteps)

# Install compiled lib
install(TARGETS cysteps LIBRARY DESTINATION "${PYTHON_INSTALL_PREFIX}/steps/")

# Install python package, from cmake_setup.py
install(
  CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${setup_py} build --build-base=${CMAKE_BINARY_DIR} install --install-lib=${PYTHON_INSTALL_PREFIX} --force WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})"
)
