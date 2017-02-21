# FindPythonInterpLibs
# --------------------
#
# Find the python libraries and include path based on
# what the python interpreter reports.
#
# In addition to the variables set by FindPythonInterp(), sets:
#    PYTHONLIBS_FOUND       Success?
#    PYTHON_LIBRARIES       Path to python library
#    PYTHON_INCLUDE_DIRS    Path to Python.h include file

if(PythonInterpLibs_FIND_REQUIRED)
    set(_required "REQUIRED")
endif()

if(PythonInterpLibs_FIND_QUIET)
    set(_quiet "QUIET")
endif()

find_package(PythonInterp ${_required} ${_quiet})

if(PYTHON_EXECUTABLE)
    set(_pylib_name "python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import distutils.sysconfig as ds; print ds.PREFIX"
                    OUTPUT_VARIABLE _prefix OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(_prefix)
        find_library(PYTHON_LIBRARIES "${_pylib_name}" PATHS "${_prefix}/lib" PATH_SUFFIXES ${CMAKE_LIBRARY_ARCHITECTURE} NO_DEFAULT_PATH)
    endif()
   
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import distutils.sysconfig as ds; print ds.get_python_inc(True)"
                    OUTPUT_VARIABLE _inc_path)
    if(_inc_path)
        set(PYTHON_INCLUDE_DIRS "${_inc_path}" CACHE PATH "Path to where Python.h is found")
    endif()
endif()
unset(_prefix)
unset(_inc_path)
unset(_pylib_name)
unset(_required)
unset(_quiet)

find_package_handle_standard_args(PythonInterpLibs DEFAULT_MSG PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS)


