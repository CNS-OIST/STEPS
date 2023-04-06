if(POLICY CMP0114)
  cmake_policy(SET CMP0114 NEW)
endif()

# set path for bundled build
set(tpdir "${CMAKE_BINARY_DIR}/_bundle")
set(Omega_h_libs "${tpdir}/omega_h-install/lib/libomega_h${CMAKE_SHARED_LIBRARY_SUFFIX}")

# reconfigure only if library is missing (otherwise it recompiles every time)
include(ExternalProject)
ExternalProject_Add(
  omega_h
  BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/steps/libomega_h${CMAKE_SHARED_LIBRARY_SUFFIX}"
  GIT_REPOSITORY https://github.com/sandialabs/omega_h.git
  GIT_TAG v9.34.13
  TMP_DIR "${tpdir}/omega_h-tmp"
  STAMP_DIR "${tpdir}/omega_h-stamp"
  DOWNLOAD_DIR "${tpdir}/omega_h-src"
  SOURCE_DIR "${tpdir}/omega_h-src"
  BINARY_DIR "${tpdir}/omega_h-build"
  INSTALL_DIR "${tpdir}/omega_h-install"
  UPDATE_DISCONNECTED TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=../omega_h-install -DOmega_h_USE_MPI:BOOL=ON
             -DCMAKE_CXX_COMPILER:FILEPATH=mpicxx -DOmega_h_USE_Gmsh_DEFAULT:BOOL=ON
             -DOmega_h_USE_Gmsh:BOOL=ON -DCMAKE_CXX_FLAGS=-O3
  EXCLUDE_FROM_ALL True)
ExternalProject_Add_StepTargets(omega_h install)

# create include path that may not exist yet
file(MAKE_DIRECTORY ${tpdir}/omega_h-install/include)

# add target
include(ImportModernLib)
add_library_target(NAME Omega_h
  INCLUDE_DIRECTORIES
    "${tpdir}/omega_h-install/include"
  LIBRARIES
    "${tpdir}/omega_h-install/lib/libomega_h${CMAKE_SHARED_LIBRARY_SUFFIX}"
  NAME_LOWER_CASE True)

# ensure that omega_h is compiled and installed before the target is linked
add_dependencies(Omega_h::omega_h omega_h-install)
