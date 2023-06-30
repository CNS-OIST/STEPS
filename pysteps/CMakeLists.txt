# ------------------------------------------------------------------
# Definitions
# ------------------------------------------------------------------
set(steps_include_dir "${PROJECT_SOURCE_DIR}/src/steps")
set(localpkgs "${CMAKE_BINARY_DIR}/packages")
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
set(steps_requirements_file "${CMAKE_BINARY_DIR}/steps-requirements.txt")

# Cython
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

# ------------------------------------------------------------------
# Make rules
# ------------------------------------------------------------------
# .pyx -> .cpp
if(PETSC_FOUND AND USE_PETSC)
  set(cython_use_petsc True)
else()
  set(cython_use_petsc False)
endif()
if(STEPS_USE_DIST_MESH)
  set(cython_use_dist_mesh True)
else()
  set(cython_use_dist_mesh False)
endif()
add_custom_command(
  OUTPUT ${cython_outputs}
  COMMAND
    ${CYTHON_EXECUTABLE} --fast-fail --cplus -I ${steps_include_dir} --compile-time-env
    USE_PETSC=${cython_use_petsc} --compile-time-env STEPS_USE_DIST_MESH=${cython_use_dist_mesh}
    --output-file ${cython_outputs} ${cython_src} -X embedsignature=True
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS libsteps_static ${cy_srcs} ${steps_package_timestamp}
  COMMENT "Cythonizing extension ${cython_src}")

# .cpp -> .so (it will be a single .so without dependencies on libsteps)
set(CYSTEPS_COMPILE_FLAGS
    "-Wno-strict-aliasing -Wno-tautological-compare \
     -Wno-format -Wno-double-promotion -Wno-unsequenced \
     -Wno-format-nonliteral \
     -Wno-unused-parameter -Wno-old-style-cast -Wno-unused-variable \
     -Wno-shadow -Wno-parentheses -Wno-cast-align")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CYSTEPS_COMPILE_FLAGS
      "${CYSTEPS_COMPILE_FLAGS} -Wno-useless-cast -Wno-deprecated-declarations -Wno-overflow")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
    set(CYSTEPS_COMPILE_FLAGS "${CYSTEPS_COMPILE_FLAGS} -Wno-cast-function-type")
  endif()
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  set(CYSTEPS_COMPILE_FLAGS
      "${CYSTEPS_COMPILE_FLAGS} -Wno-c++17-extensions -Wno-implicit-fallthrough")
endif()

# Set RPATH options
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${bundle_lib_install_dir}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(APPLE)
  # so that the cython library can find the bundled shared libraries copied there.
  list(APPEND CMAKE_BUILD_RPATH ${CMAKE_BINARY_DIR}/lib/steps)
endif(APPLE)

python_add_module(cysteps ${cython_outputs})

set_source_files_properties(${cython_outputs} PROPERTIES COMPILE_FLAGS "${CYSTEPS_COMPILE_FLAGS}")

target_include_directories(cysteps PRIVATE ${PYTHON_INCLUDE_DIRS} ${steps_include_dir})
target_link_libraries(cysteps libsteps_static ${libsteps_link_libraries})

set_target_properties(cysteps PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${localpkgs}/steps" OUTPUT_NAME
                                                                                       ${cymodname})

set_directory_properties(
  PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
             "${CMAKE_CURRENT_SOURCE_DIR}/build;${CMAKE_CURRENT_SOURCE_DIR}/steps.egg-info")

# Python optional packages and extras
set(STEPS_INCLUDE_PACKAGES "steps" "steps.*")
if(STEPS_USE_DIST_MESH)
  list(APPEND py_deps "dist")
endif()
if(STEPS_USE_HDF5_SAVING)
  list(APPEND py_deps "hdf5")
endif()
if(ENABLE_PYTHON_CODECOVERAGE)
  list(APPEND py_deps "cov")
endif()
if(BUILD_TESTING)
  list(APPEND py_deps "test")
endif()
if(STEPS_USE_STEPSBLENDER)
  list(APPEND py_deps "blender")
  list(APPEND STEPS_INCLUDE_PACKAGES "stepsblender*")
endif()
list(JOIN py_deps "," py_all_deps)
list(JOIN STEPS_INCLUDE_PACKAGES "," STEPS_INCLUDE_PACKAGES_SETUP)
list(JOIN STEPS_INCLUDE_PACKAGES "\",\"" STEPS_INCLUDE_PACKAGES_PYPROJECT)
set(STEPS_INCLUDE_PACKAGES_PYPROJECT "\"${STEPS_INCLUDE_PACKAGES_PYPROJECT}\"")

# Configure pyproject.toml and setup.cfg to only include required subpackages
configure_file(pyproject.toml.in "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml")
configure_file(setup.cfg.in "${CMAKE_CURRENT_SOURCE_DIR}/setup.cfg")

add_custom_command(
  OUTPUT ${steps_package_timestamp}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND ${PYTHON_EXECUTABLE} -m build --wheel --outdir ${CMAKE_BINARY_DIR}
  COMMAND ${PYTHON_EXECUTABLE} -m pip install -t ${localpkgs} --upgrade --no-deps --find-links
          ${CMAKE_BINARY_DIR} --no-index steps[${py_all_deps}]
  COMMAND ${CMAKE_COMMAND} -E touch ${steps_package_timestamp})

# Install Python package with pip
if(STEPS_INSTALL_PYTHON_DEPS)
  # Pure python package
  add_custom_target(
    steps_package_pure_py ALL DEPENDS ${steps_package_timestamp} ${steps_requirements_file}
  )# random file bridging
  # target->command

  if(PIPTOOLS_VERSION_STRING VERSION_LESS 6.8.0)
    add_custom_command(
      OUTPUT ${steps_requirements_file}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${PYTHON_EXECUTABLE} -m piptools compile --quiet --extra ${py_all_deps} -o
              ${steps_requirements_file} pyproject.toml)
  else()
    add_custom_command(
      OUTPUT ${steps_requirements_file}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${PYTHON_EXECUTABLE} -m piptools compile --quiet --extra ${py_all_deps} -o
              ${steps_requirements_file} --resolver=backtracking pyproject.toml)
  endif()

  install(
    CODE "execute_process(
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${PYTHON_EXECUTABLE} -m build --wheel --outdir ${CMAKE_BINARY_DIR})
  execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -m pip install --requirement ${steps_requirements_file})
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m pip install -t ${PYTHON_INSTALL_PREFIX} --upgrade --force-reinstall --no-deps --find-links ${CMAKE_BINARY_DIR} --no-index steps[${py_all_deps}])"
  )
else()
  install(
    CODE "execute_process(
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${PYTHON_EXECUTABLE} -m build --wheel --outdir ${CMAKE_BINARY_DIR})
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m pip install -t ${PYTHON_INSTALL_PREFIX} --upgrade --force-reinstall --no-deps --find-links ${CMAKE_BINARY_DIR} --no-index steps[${py_all_deps}])"
  )

endif()

# Install compiled lib
install(TARGETS cysteps LIBRARY DESTINATION "${PYTHON_INSTALL_PREFIX}/steps/")
