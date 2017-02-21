file(COPY steps DESTINATION .)

# For 'build' step, install steps package under build/lib directory.
# For 'install' step, perform distutils install under CMAKE_INSTALL_PREFIX

set(lib_dir "${CMAKE_BINARY_DIR}/lib")

file(GLOB_RECURSE steps_package_sources steps *.py)

set(steps_package_timestamp "${CMAKE_CURRENT_BINARY_DIR}/steps_package_timestamp")
set(steps_packages "'steps','steps/utilities','steps/visual'")
if(MPI_FOUND)
    set(steps_packages "${steps_packages},'steps/mpi'")
endif()
configure_file(cmake_setup.py.in cmake_setup.py)

# Variables steps_swig_target and steps_swig_pymodules are set in swig/CMakeLists.txt,
# and may depend on the availability of numpy.

add_custom_target(steps_package_rebuild_check ALL DEPENDS ${steps_package_timestamp})
add_custom_command(OUTPUT ${steps_package_timestamp}
    COMMAND ${CMAKE_COMMAND} -E copy ${steps_swig_pymodules} steps
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${steps_swig_target}> steps
    COMMAND ${PYTHON_EXECUTABLE} cmake_setup.py install_lib "--install-dir" "${lib_dir}" "--force"
    COMMAND ${CMAKE_COMMAND} -E touch ${steps_package_timestamp}
    DEPENDS ${steps_swig_target} ${steps_package_sources})

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} cmake_setup.py install --prefix ${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")

