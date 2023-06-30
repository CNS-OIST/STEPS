# Find Python numpy module.
# ~~~
# Will set:
# NUMPY_FOUND if found;
# NUMPY_VERSION_STRING if numpy module __version__ is defined;
# NUMPY_INCLUDE_DIRS to the include directory reported by numpy.
# ~~~

include(FindPackageHandleStandardArgs)

if(NumPy_FIND_QUIETLY)
  set(_NUMPY_QUIET "QUIET")
else()
  unset(_NUMPY_QUIET)
endif()

if(NumPy_FIND_REQUIRED)
  set(_NUMPY_REQUIRED "REQUIRED")
else()
  unset(_NUMPY_REQUIRED)
endif()

find_package(PythonInterp ${_NUMPY_QUIET} ${_NUMPY_REQUIRED})

if(NOT PYTHON_EXECUTABLE)
  fail_message("No python interpreter found")
  set(NUMPY_FOUND 0)
else()
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
            "import numpy; print numpy.__version__ if hasattr(numpy,'__version__') else ''"
    ERROR_QUIET
    RESULT_VARIABLE _NUMPY_RV
    OUTPUT_VARIABLE NUMPY_VERSION_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_NUMPY_RV)
    unset(NUMPY_VERSION_STRING)
  endif()

  find_package_handle_standard_args(
    NumPy
    REQUIRED_VARS NUMPY_VERSION_STRING
    VERSION_VAR NUMPY_VERSION_STRING)

  if(NUMPY_FOUND)
    set(_NUMPY_GET_INCLUDE "get_include")
    if(NUMPY_VERSION_STRING VERSION_LESS "1.0")
      set(_NUMPY_GET_INCLUDE "get_numpy_include")
    endif()

    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.${_NUMPY_GET_INCLUDE}()"
      ERROR_QUIET
      OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS)
  endif()
endif()
