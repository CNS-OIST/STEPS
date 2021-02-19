# Patch to use openmp together with AppleClang and old cmake style.
#
# It is not needed for modern cmake (>=3.12), it does not add the library based
# on target.
#
# It tries to find the OpenMP library installed by brew, anaconda, or macports
#

message(STATUS "Applying patch for OpenMP and AppleClang")

# Try brew
find_program(BREW_EXECUTABLE "brew")
mark_as_advanced(BREW_EXECUTABLE)
if(BREW_EXECUTABLE)
  execute_process(COMMAND bash "-c" "brew --prefix libomp" RESULT_VARIABLE ACP_OMP_LIB_FOUND OUTPUT_VARIABLE ACP_OMP_BASE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(ACP_OMP_LIB_FOUND EQUAL 0)
    set(ACP_OMP_INCLUDE_PATH ${ACP_OMP_BASE_PATH}/include)
    set(ACP_OMP_LIB_PATH ${ACP_OMP_BASE_PATH}/lib)
    set(ACP_OMP_LIB ${ACP_OMP_LIB_PATH}/libomp.dylib)
  endif()
endif()

# Try anaconda only if brew failed
find_program(ANACONDA_EXECUTABLE "conda")
mark_as_advanced(ANACONDA_EXECUTABLE)
if(ANACONDA_EXECUTABLE AND NOT EXISTS ${ACP_OMP_LIB})
  execute_process(COMMAND bash "-c" "echo $CONDA_PREFIX" RESULT_VARIABLE CONDA_PREFIX_FOUND OUTPUT_VARIABLE ACP_OMP_BASE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(CONDA_PREFIX_FOUND EQUAL 0)
    set(ACP_OMP_INCLUDE_PATH ${ACP_OMP_BASE_PATH}/include)
    set(ACP_OMP_LIB_PATH ${ACP_OMP_BASE_PATH}/lib)
    set(ACP_OMP_LIB ${ACP_OMP_LIB_PATH}/libomp.dylib)
  endif()
endif()

# Try macports only if anaconda and brew failed
find_program(MACPORT_EXECUTABLE "port")
mark_as_advanced(MACPORT_EXECUTABLE)
if(MACPORT_EXECUTABLE AND NOT EXISTS ${ACP_OMP_LIB})
  execute_process(COMMAND bash "-c" "port -q contents libomp | grep libomp.dylib" RESULT_VARIABLE ACP_OMP_LIB_FOUND OUTPUT_VARIABLE ACP_OMP_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(ACP_OMP_LIB_FOUND EQUAL 0)
    # Strip leading spaces
    string(STRIP ${ACP_OMP_LIB} ACP_OMP_LIB)
    # Strip /libomp.dylib
    get_filename_component(ACP_OMP_LIB_PATH ${ACP_OMP_LIB} DIRECTORY)
    # Strip /libomp
    get_filename_component(ACP_OMP_BASE_PATH ${ACP_OMP_LIB_PATH} DIRECTORY)
    # Strip /lib
    get_filename_component(ACP_OMP_BASE_PATH ${ACP_OMP_BASE_PATH} DIRECTORY)
    set(ACP_OMP_INCLUDE_PATH ${ACP_OMP_BASE_PATH}/include/libomp)
  endif()
endif()

# Set flags, include dirs, and lib
if(EXISTS ${ACP_OMP_LIB})
  message(STATUS "Found libomp: ${ACP_OMP_LIB}")
  set(OPENMP_FOUND True)

  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "CXX compiler flags for OpenMP parallelization" FORCE)
  set(OpenMP_CXX_INCLUDE_DIR ${ACP_OMP_INCLUDE_PATH} CACHE STRING "Include path" FORCE)
  set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "CXX compiler libraries for OpenMP parallelization" FORCE)

  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "C compiler flags for OpenMP parallelization" FORCE)
  set(OpenMP_C_INCLUDE_DIR ${ACP_OMP_INCLUDE_PATH} CACHE STRING "Include path" FORCE)
  set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "C compiler libraries for OpenMP parallelization" FORCE)

  set(OpenMP_libomp_LIBRARY ${ACP_OMP_LIB} CACHE STRING "Path to library." FORCE)

else()
  set(OPENMP_FOUND False)
  message(WARNING "Could NOT find libomp ")
endif()