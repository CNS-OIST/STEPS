if(STEPS_USE_DIST_MESH)
    add_test(NAME python-unit-tests-dist
            COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
            ${PYTHON_EXECUTABLE} run.py
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME python-unit-tests-dist-2
            COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${PROJECT_BINARY_DIR}/lib:$ENV{PYTHONPATH}"
            ${MPIRUN} -n 2 ${PYTHON_EXECUTABLE} run.py
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
