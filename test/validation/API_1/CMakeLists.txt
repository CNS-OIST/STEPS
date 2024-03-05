add_python_tests(python-validation-API_1 . ${SERIAL_PYTEST_PREFIX})

if(MPI_FOUND)
  add_python_tests(python-validation-parallel-API_1 . ${PARALLEL_PYTEST_PREFIX})
endif()
