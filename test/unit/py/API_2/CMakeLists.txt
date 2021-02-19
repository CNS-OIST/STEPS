add_test(NAME python-unit-tests-API2
         COMMAND ${CMAKE_COMMAND} -E env
                 "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
                 ${PYTHON_EXECUTABLE} run_unit_tests.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
if(MPI_FOUND)
    add_test(NAME parallel-python-unit-tests-API2
           COMMAND ${CMAKE_COMMAND} -E env
                   "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
                   ${MPIRUN} -n 2 ${PYTHON_EXECUTABLE}
                   run_parallel_unit_tests.py
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
