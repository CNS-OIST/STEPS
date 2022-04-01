# - Add python tests that use the 'unittest' module
#
# Synopsis:
#
#	add_python_tests (path prefix)
#
# where:
#
#   name       Name of the group of tests, it will be used as part of the names of individual
#              tests. In general the name of tests will be of the form: 
#              "{name}.{testCase}.{testMethod}" with {testCase} the name of the class inheriting
#              from unittest.TestCase and {testMethod} the name of the test method in this class.
#
#   path       Path to a directory in which the tests should be found. This value is
#              passed to unittest's discover method as 'start_dir'.
#
#   prefix     Only tests located in files whose name match '${prefix}*.py'
#              will be added. Only ${SERIAL_PYTEST_PREFIX}, ${PARALLEL_PYTEST_PREFIX} and
#              ${DISTRIBUTED_PYTEST_PREFIX} are accepted values.

set(SERIAL_PYTEST_PREFIX test_)
set(PARALLEL_PYTEST_PREFIX parallel_)
set(DISTRIBUTED_PYTEST_PREFIX distributed_)

function(add_python_tests name path prefix)
    # List tests using unittest's 'discover' method
    execute_process(COMMAND ${CMAKE_COMMAND} -E env
                            "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
                            ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/CMake/list_tests.py "${path}" "${prefix}"
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    OUTPUT_VARIABLE STR_TESTS_PATHS
                    ERROR_VARIABLE STR_LOADING_ERRORS
                    RESULT_VARIABLE EXIT_CODE
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(EXIT_CODE EQUAL "1")
        message(WARNING "While loading ${name} tests: ${STR_LOADING_ERRORS}")
    endif()

    # Add all discovered tests
    foreach(TEST_PATH ${STR_TESTS_PATHS})
        string(REGEX MATCH "[A-Za-z0-9_]+\\.[A-Za-z0-9_]+$" TEST_NAME "${TEST_PATH}")
        if(${prefix} STREQUAL ${SERIAL_PYTEST_PREFIX})
            add_test(NAME "${name}.${TEST_NAME}"
                     COMMAND ${CMAKE_COMMAND} -E env
                             "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
                             ${PYTHON_EXECUTABLE} -m unittest ${TEST_PATH}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        elseif(${prefix} STREQUAL ${PARALLEL_PYTEST_PREFIX})
            add_test(NAME "${name}.${TEST_NAME}"
                     COMMAND ${CMAKE_COMMAND} -E env
                             "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
                             ${MPIRUN} -n 2 ${PYTHON_EXECUTABLE} -m unittest ${TEST_PATH}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        elseif(${prefix} STREQUAL ${DISTRIBUTED_PYTEST_PREFIX})
            add_test(NAME "${name}.${TEST_NAME}"
                     COMMAND ${CMAKE_COMMAND} -E env
                             "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
                             ${MPIRUN} -n 2 ${PYTHON_EXECUTABLE} -m unittest ${TEST_PATH}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        else()
            message(WARNING "Incorrect python test prefix used, corresponding tests were not added.")
        endif()
    endforeach(TEST_PATH)
endfunction()
