add_python_tests(python-unit-API_1 . ${SERIAL_PYTEST_PREFIX})

if(MPI_FOUND)
    add_python_tests(python-unit-parallel-API_1 . ${PARALLEL_PYTEST_PREFIX})
endif()
