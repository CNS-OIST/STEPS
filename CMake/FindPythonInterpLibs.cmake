# Finds Python libs and paths
# ---------------------------
#  This module finds the Site-packages, Include dir and shared lib of the Python interpreter.
#  It also provides PYTHON_ADD_MODULE function to define Python extension modules
#
#  PYTHON_SITE_PACKAGES    = path to the python modules dir
#  PYTHON_LIBRARIES        = path to the python library
#  PYTHON_INCLUDE_DIRS     = path to where Python.h is found
# --
find_package(PythonInterp REQUIRED)

## Python PREFIX
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import EXEC_PREFIX; print(EXEC_PREFIX)"
                OUTPUT_VARIABLE PYTHON_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)

## Site-packages
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(True))"
                OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
                OUTPUT_STRIP_TRAILING_WHITESPACE)

## Include dir
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"
                OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

## Lib and devel dirs
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; \
                        print('%s;%s' % (sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LIBPL')))"
                OUTPUT_VARIABLE PYTHON_LIBDIRS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

## Find lib name
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; \
                        print('%s;%s' % (sysconfig.get_config_var('LDLIBRARY'), sysconfig.get_config_var('LIBRARY')))"
                OUTPUT_VARIABLE PYTHON_LIB_NAMES
                OUTPUT_STRIP_TRAILING_WHITESPACE)
                

find_library(PYTHON_LIBRARY
    #NAMES ${_lib_name} ${_lib_name}mu ${_lib_name}m ${_lib_name}u
    NAMES ${PYTHON_LIB_NAMES}
    PATHS ${PYTHON_LIBDIRS}
    NO_DEFAULT_PATH)
    
if(PYTHON_LIBRARY MATCHES ".a$")
    message(WARNING "Python only found as a static library")
endif()


## Setting cached vars
set(PYTHON_SITE_PACKAGES ${PYTHON_SITE_PACKAGES}
    CACHE PATH "path to the python modules dir")

set(PYTHON_LIBRARIES ${PYTHON_LIBRARY}
    CACHE PATH "path to the python library")

set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS}
    CACHE PATH "path to the python include dir")

mark_as_advanced(PYTHON_SITE_PACKAGES)
mark_as_advanced(PYTHON_LIBRARIES)
mark_as_advanced(PYTHON_INCLUDE_DIRS)

# Handle Options
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PythonInterpLibs DEFAULT_MSG PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS)


# Make default destination for Python packages
set(_default_prefixes "/usr/local" "/usr")
list(FIND _default_prefixes ${CMAKE_INSTALL_PREFIX} _index)
if(${_index} GREATER -1)
    # Default installations install Python to default site-packages
    set(_py_install_prefix ${PYTHON_SITE_PACKAGES})
else()
    set(_py_install_prefix ${CMAKE_INSTALL_PREFIX})
endif()

set(PYTHON_INSTALL_PREFIX ${_py_install_prefix} CACHE PATH "Target install prefix for Python modules")
mark_as_advanced(PYTHON_INSTALL_PREFIX)


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
    # Clang requires explicitly setting dynamic_lookup on undefined symbols
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_link_libraries(${_NAME} "-Wl,-undefined,dynamic_lookup")
    endif()

    if(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      set_target_properties(${_NAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
      if(WIN32 AND NOT CYGWIN)
        set_target_properties(${_NAME} PROPERTIES SUFFIX ".pyd")
      endif()
    endif()

  endif()
endfunction()
