if(POLICY CMP0114)
  cmake_policy(SET CMP0114 NEW)
endif()

# set path for bundled build
set(tpdir "${CMAKE_BINARY_DIR}/_bundle")
set(SUNDIALS_libs "${tpdir}/SUNDIALS-install/lib/libsundials_cvode${CMAKE_SHARED_LIBRARY_SUFFIX}")

# reconfigure only if library is missing (otherwise it recompiles every time)
include(ExternalProject)
ExternalProject_Add(
  SUNDIALS
  BUILD_BYPRODUCTS
    "${CMAKE_BINARY_DIR}/lib/steps/libsundials_nvecserial${CMAKE_SHARED_LIBRARY_SUFFIX}"
    "${CMAKE_BINARY_DIR}/lib/steps/libsundials_cvode${CMAKE_SHARED_LIBRARY_SUFFIX}"
    "${CMAKE_BINARY_DIR}/lib/steps/libsundials_ida${CMAKE_SHARED_LIBRARY_SUFFIX}"
    "${CMAKE_BINARY_DIR}/lib/steps/libsundials_kinsol${CMAKE_SHARED_LIBRARY_SUFFIX}"
  DOWNLOAD_NO_PROGRESS TRUE
  URL "https://github.com/LLNL/sundials/releases/download/v3.2.1/sundials-3.2.1.tar.gz"
  URL_HASH SHA256=47d94d977ab2382cdcdd02f72a25ebd4ba8ca2634bbb2f191fe1636e71c86808
  TMP_DIR "${tpdir}/SUNDIALS-tmp"
  STAMP_DIR "${tpdir}/SUNDIALS-stamp"
  DOWNLOAD_DIR "${tpdir}/SUNDIALS-src"
  SOURCE_DIR "${tpdir}/SUNDIALS-src"
  BINARY_DIR "${tpdir}/SUNDIALS-build"
  INSTALL_DIR "${tpdir}/SUNDIALS-install"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=../SUNDIALS-install -DBUILD_IDAS=OFF -DBUILD_CVODES=OFF
             -DBUILD_STATIC_LIBS=OFF -DCMAKE_INSTALL_LIBDIR=lib
  EXCLUDE_FROM_ALL True)
ExternalProject_Add_StepTargets(SUNDIALS install)

# create include path that may not exist yet
file(MAKE_DIRECTORY ${tpdir}/SUNDIALS-install/include)

# add target
include(ImportModernLib)
add_library_target(NAME SUNDIALS
  INCLUDE_DIRECTORIES "${tpdir}/SUNDIALS-install/include"
  LIBRARIES "${tpdir}/SUNDIALS-install/lib/libsundials_nvecserial${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${tpdir}/SUNDIALS-install/lib/libsundials_cvode${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${tpdir}/SUNDIALS-install/lib/libsundials_ida${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${tpdir}/SUNDIALS-install/lib/libsundials_kinsol${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

# ensure that SUNDIALS is compiled and installed before the target is linked
add_dependencies(SUNDIALS::SUNDIALS SUNDIALS-install)

if(CMAKE_VERSION VERSION_LESS 3.12.0)
  add_definitions(-DSTEPS_SUNDIALS_VERSION_MAJOR=3)
else()
  add_compile_definitions(STEPS_SUNDIALS_VERSION_MAJOR=3)
endif()
