# - Find SWIG

find_program(SWIG_EXECUTABLE NAMES swig3.0 swig2.0 swig)

if(SWIG_EXECUTABLE)
    execute_process(COMMAND ${SWIG_EXECUTABLE} -swiglib
        OUTPUT_VARIABLE SWIG_DIR RESULT_VARIABLE _RV ERROR_QUIET)
    if(_RV)
	set(SWIG_DIR SWIG_DIR-NOTFOUND)
    endif()

    if(SWIG_DIR)
        execute_process(COMMAND ${SWIG_EXECUTABLE} -version
            OUTPUT_VARIABLE SWIG_VERSION RESULT_VARIABLE _RV ERROR_QUIET)

        if(_RV)
            unset(SWIG_VERSION)
        else()
            string(REGEX REPLACE ".*SWIG Version ([0-9.]+).*" "\\1" SWIG_VERSION "${SWIG_VERSION}")
        endif()
    endif()
endif()

unset(_RV)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SWIG REQUIRED_VARS SWIG_EXECUTABLE SWIG_DIR VERSION_VAR SWIG_VERSION)
