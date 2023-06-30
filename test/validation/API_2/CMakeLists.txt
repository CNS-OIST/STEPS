add_python_tests(python-validation-API_2 . ${SERIAL_PYTEST_PREFIX})

if(MPI_FOUND)
  add_python_tests(python-validation-parallel-API_2 . ${PARALLEL_PYTEST_PREFIX})
endif()
