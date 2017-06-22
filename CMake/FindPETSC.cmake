# Try to find PETSC, without performing long tests
# System directories and additional CMAKE_PREFIX_PATH are already searched automatically
find_library(PETSC_LIBRARIES
  NAMES petsc
  HINTS "$ENV{PETSC_DIR}/lib"
)
find_path(PETSC_INCLUDE_DIRS
  NAMES petsc.h
  HINTS "$ENV{PETSC_DIR}/include"
  PATH_SUFFIXES petsc
)

# Checks argurments 'REQUIRED', 'QUIET', and version.
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(PETSC FOUND_VAR PETSC_FOUND REQUIRED_VARS PETSC_LIBRARIES PETSC_INCLUDE_DIRS)
