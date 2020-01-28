# ------------------------------------------------------------------
# Definitions
# ------------------------------------------------------------------
set(steps_include_dir "${PROJECT_SOURCE_DIR}/src")
set(setup_py "${CMAKE_CURRENT_BINARY_DIR}/cmake_setup.py")
set(libdir "${CMAKE_BINARY_DIR}/lib")
set(cymodname cysteps)
if(MPI_FOUND)
  set(cymodname ${cymodname}_mpi)
endif()
set(cython_output "${CMAKE_CURRENT_SOURCE_DIR}/steps/${cymodname}.cpp")
set(steps_package_timestamp "${CMAKE_CURRENT_BINARY_DIR}/python_built.tmp")
set(cython_src ${cymodname}.pyx)
# Track cython sources
file(GLOB cy_srcs *.pyx *.pxd)

# Dont add build to RPATH
set(CMAKE_SKIP_BUILD_RPATH TRUE)

# Configure setup.py to include version
configure_file(cmake_setup.py.in ${setup_py})

# ------------------------------------------------------------------
# Make rules
# ------------------------------------------------------------------
# .pyx -> .cpp
add_custom_command(OUTPUT ${cython_output}
                   COMMAND ${CYTHON_EXECUTABLE}
                           --fast-fail
                           --cplus -I ${steps_include_dir}
                           --output-file ${cython_output} ${cython_src}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   DEPENDS libsteps ${cy_srcs}
                   COMMENT "Cythonizing extension ${cython_src}")


# .cpp -> .so (it will be a single .so without dxependencies on libsteps)
set(CYSTEPS_COMPILE_FLAGS
    "-Wno-strict-aliasing -Wno-tautological-compare -Wno-useless-cast \
     -Wno-format -Wno-double-promotion -Wno-unsequenced \
     -Wno-unused-parameter -Wno-old-style-cast -Wno-unused-variable \
     -Wno-shadow -Wno-parentheses -Wno-cast-align")

python_add_module(cysteps ${cython_output})

set_source_files_properties(${cython_output}
        PROPERTIES
        COMPILE_FLAGS "${CYSTEPS_COMPILE_FLAGS}")

target_include_directories(cysteps
                           PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI_FOUND_HEADERS}
                                   ${steps_include_dir})
target_link_libraries(cysteps libsteps_static ${libsteps_link_libraries})

set_target_properties(cysteps
                      PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                 "${libdir}/steps"
                                 OUTPUT_NAME
                                 ${cymodname})

# Pure python package
add_custom_target(steps_package_pure_py ALL
                  DEPENDS ${steps_package_timestamp}) # random file bridging
                                                      # target->command
add_custom_command(OUTPUT ${steps_package_timestamp}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   COMMAND ${PYTHON_EXECUTABLE} ${setup_py} build
                           --build-base=${CMAKE_BINARY_DIR} install
                           --install-lib "${libdir}"
                           --no-compile
                           --force
                   COMMAND ${CMAKE_COMMAND} -E touch ${steps_package_timestamp}
                           # setup.py already tracks changed .py files. We set
                           # DEPENDS just to set the order: copy .py after
                           # building the lib
                   DEPENDS cysteps)

# Install compiled lib
install(TARGETS cysteps LIBRARY DESTINATION "${PYTHON_INSTALL_PREFIX}/steps/")

# Install python package, from cmake_setup.py
install(
  CODE
  "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${setup_py} build --build-base=${CMAKE_BINARY_DIR} install --install-lib=${PYTHON_INSTALL_PREFIX} --force WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})"
  )
