# .rst: FindGmsh
# ~~~
# --------------
#
# Find Gmsh executable, headers and library
#
# ::
#
# * Gmsh_FOUND           - True if Gmsh library is found.
# * Gmsh_INCLUDE_DIRS    - Directory where Gmsh headers are located.
# * Gmsh_LIBRARIES       - Gmsh libraries to link against.
# * Gmsh_VERSION_MAJOR   - The major version of Gmsh.
# * Gmsh_VERSION_MINOR   - The minor version of Gmsh.
# * Gmsh_VERSION_PATCH   - The patch version of Gmsh.
# * Gmsh_VERSION_STRING  - version number as a string (ex: "4.3.2").
# ~~~

find_program(Gmsh_EXECUTABLE gmsh)

# Extract version from command "gmsh --version"
if(Gmsh_EXECUTABLE)
  execute_process(
    COMMAND ${Gmsh_EXECUTABLE} --version
    ERROR_VARIABLE Gmsh_VERSION_STRING
    OUTPUT_VARIABLE Gmsh_VERSION_STRING
    ERROR_STRIP_TRAILING_WHITESPACE)

  if(Gmsh_VERSION_STRING)
    # Gmsh_VERSION sample: "4.4.1-git-ba5563507"
    string(REGEX REPLACE "^([.0-9]+).*" "\\1" Gmsh_VERSION_STRING "${Gmsh_VERSION_STRING}")
    # Extract version components
    string(REPLACE "." ";" Gmsh_version "${Gmsh_VERSION_STRING}")
    list(LENGTH Gmsh_version Gmsh_VERSION_COUNT)
    if(Gmsh_VERSION_COUNT GREATER 0)
      list(GET Gmsh_version 0 Gmsh_VERSION_MAJOR)
    else()
      set(Gmsh_VERSION_MAJOR 0)
    endif()
    if(Gmsh_VERSION_COUNT GREATER 1)
      list(GET Gmsh_version 1 Gmsh_VERSION_MINOR)
    else()
      set(Gmsh_VERSION_MINOR 0)
    endif()
    if(Gmsh_VERSION_COUNT GREATER 2)
      list(GET Gmsh_version 2 Gmsh_VERSION_PATCH)
    else()
      set(Gmsh_VERSION_PATCH 0)
    endif()
    unset(Gmsh_version)
  endif()
endif()

# First try to find include and libraries from pip-installed gmsh
execute_process(
  COMMAND "${Python_EXECUTABLE}" "-c" "import gmsh; print(gmsh.libpath, end='')"
  ERROR_VARIABLE Gmsh_PIP_LIBRARIES
  OUTPUT_VARIABLE Gmsh_PIP_LIBRARIES
  ERROR_STRIP_TRAILING_WHITESPACE)
if(EXISTS "${Gmsh_PIP_LIBRARIES}")
  set(Gmsh_LIBRARIES "${Gmsh_PIP_LIBRARIES}")
  cmake_path(GET Gmsh_PIP_LIBRARIES PARENT_PATH Gmsh_PIP_LIBDIR)
  find_path(
    Gmsh_INCLUDE_DIRS
    NAMES gmsh.h
    HINTS "${Gmsh_PIP_LIBDIR}/../include")
else()
  # If gmsh internals changed and gmsh.libpath is no longer accessible, try to find the includes and
  # library with find_path and find_library

  # First try to find a local gmsh (e.g. installed from pip)
  find_path(Gmsh_INCLUDE_DIRS NAMES gmsh.h NO_CMAKE_SYSTEM_PATH)
  find_library(
    Gmsh_LIBRARIES
    NAMES gmsh gmsh.${Gmsh_VERSION_MAJOR}.${Gmsh_VERSION_MINOR}
          libgmsh.so.${Gmsh_VERSION_MAJOR}.${Gmsh_VERSION_MINOR} NO_CMAKE_SYSTEM_PATH)

  # If that didn't work, look for system files as well
  if(Gmsh_INCLUDE_DIRS STREQUAL "Gmsh_INCLUDE_DIRS-NOTFOUND")
    unset(Gmsh_INCLUDE_DIRS)
    unset(Gmsh_LIBRARIES)
    find_path(Gmsh_INCLUDE_DIRS NAMES gmsh.h)
    find_library(Gmsh_LIBRARIES NAMES gmsh gmsh.${Gmsh_VERSION_MAJOR}.${Gmsh_VERSION_MINOR}
                                      libgmsh.so.${Gmsh_VERSION_MAJOR}.${Gmsh_VERSION_MINOR})
  elseif(Gmsh_LIBRARIES STREQUAL "Gmsh_LIBRARIES-NOTFOUND")
    unset(Gmsh_INCLUDE_DIRS)
    unset(Gmsh_LIBRARIES)
    find_path(Gmsh_INCLUDE_DIRS NAMES gmsh.h)
    find_library(Gmsh_LIBRARIES NAMES gmsh gmsh.${Gmsh_VERSION_MAJOR}.${Gmsh_VERSION_MINOR}
                                      libgmsh.so.${Gmsh_VERSION_MAJOR}.${Gmsh_VERSION_MINOR})
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Gmsh
  FOUND_VAR Gmsh_FOUND
  REQUIRED_VARS Gmsh_LIBRARIES Gmsh_INCLUDE_DIRS
  VERSION_VAR Gmsh_VERSION_STRING)

mark_as_advanced(Gmsh_EXECUTABLE Gmsh_INCLUDE_DIR Gmsh_LIBRARIES)
