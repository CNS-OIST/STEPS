add_test(
  NAME python-validation-tests-API2
  COMMAND
    bash -c
    "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}  ${PYTHON_EXECUTABLE} run_validation_tests.py"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
