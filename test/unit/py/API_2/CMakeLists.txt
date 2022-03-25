add_python_tests(python-unit-API_2 . ${SERIAL_PYTEST_PREFIX})

if(MPI_FOUND)
    add_python_tests(python-unit-parallel-API_2 . ${PARALLEL_PYTEST_PREFIX})
endif()

if(STEPS_USE_DIST_MESH)
    add_python_tests(python-unit-dist-API_2 . ${DISTRIBUTED_PYTEST_PREFIX})
endif()
