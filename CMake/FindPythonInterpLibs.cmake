# - Find python libraries
#  This module find the current version of python on your installation in a easy way
#
#  PYTHON_EXECUTABLE 				= python interpreter
#  PYTHON_VERSION                   = version number (MAJOR.MINOR e.g.: 2.7)
#  PYTHON_SITE_PACKAGES             = path to the python modules dir
#  PYTHON_LIBRARIES                 = path to the python library
#  PYTHON_INCLUDE_DIRS              = path to where Python.h is found
# --

include(FindPackageHandleStandardArgs)

find_package(PythonInterp REQUIRED)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%s.%s' % sys.version_info[:2])"
                                    OUTPUT_VARIABLE PYTHON_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(True))"
                OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# find include
execute_process(COMMAND
    ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"
    OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Find lib compatible with virtualenv (PYTHON_LIBRARY shall not be used in linking extensions, since they should link at runtime with the iterpreter!)
set(_pylib_name "python${PYTHON_VERSION}")
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_config_var; print(get_config_var('LIBPL'))"
                OUTPUT_VARIABLE _LIBPL OUTPUT_STRIP_TRAILING_WHITESPACE)
find_library(PYTHON_LIBRARY "${_pylib_name}" PATHS ${_LIBPL} PATH_SUFFIXES ${CMAKE_LIBRARY_ARCHITECTURE} NO_DEFAULT_PATH )

unset(_pylib_name)


# Setting cached vars
set(PYTHON_SITE_PACKAGES ${PYTHON_SITE_PACKAGES}
        CACHE PATH "path to the python modules dir")

set(PYTHON_LIBRARIES ${PYTHON_LIBRARY}
        CACHE PATH "path to the python library")

set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS}
        CACHE PATH "path to the python include dir")

mark_as_advanced(PYTHON_SITE_PACKAGES)
mark_as_advanced(PYTHON_LIBRARIES)
mark_as_advanced(PYTHON_INCLUDE_DIRS)

# Handle error
find_package_handle_standard_args(PythonLibDev DEFAULT_MSG PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS)


# function taken from the official FindPythonLibs.cmake (https://github.com/Kitware)
# PYTHON_ADD_MODULE(<name> src1 src2 ... srcN) is used to build modules for python.
# PYTHON_WRITE_MODULES_HEADER(<filename>) writes a header file you can include
# in your sources to initialize the static python modules
function(PYTHON_ADD_MODULE _NAME )
  get_property(_TARGET_SUPPORTS_SHARED_LIBS
    GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
  option(PYTHON_ENABLE_MODULE_${_NAME} "Add module ${_NAME}" TRUE)
  option(PYTHON_MODULE_${_NAME}_BUILD_SHARED "Add module ${_NAME} shared" ${_TARGET_SUPPORTS_SHARED_LIBS})

  # Mark these options as advanced
  mark_as_advanced(
    PYTHON_ENABLE_MODULE_${_NAME}
    PYTHON_MODULE_${_NAME}_BUILD_SHARED)

  if(PYTHON_ENABLE_MODULE_${_NAME})
    if(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      set(PY_MODULE_TYPE MODULE)
    else()
      set(PY_MODULE_TYPE STATIC)
      set_property(GLOBAL  APPEND  PROPERTY  PY_STATIC_MODULES_LIST ${_NAME})
    endif()

    set_property(GLOBAL  APPEND  PROPERTY  PY_MODULES_LIST ${_NAME})
    add_library(${_NAME} ${PY_MODULE_TYPE} ${ARGN})
#    target_link_libraries(${_NAME} ${PYTHON_LIBRARIES})

    if(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      set_target_properties(${_NAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
      if(WIN32 AND NOT CYGWIN)
        set_target_properties(${_NAME} PROPERTIES SUFFIX ".pyd")
      endif()
    endif()

  endif()
endfunction()
