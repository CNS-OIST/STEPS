if(POLICY CMP0114)
  cmake_policy(SET CMP0114 NEW)
endif()

# set path for bundled build
set(tpdir "${CMAKE_BINARY_DIR}/_bundle")

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
  URL "https://github.com/LLNL/sundials/releases/download/v6.0.0/sundials-6.0.0.tar.gz"
  URL_HASH SHA256=c7178e54df20a9363ae3e5ac5b3ee9db756a4ddd4b8fff045127e93b73b151f4
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

foreach(libname cvode ida kinsol nvecserial)
  if(APPLE)
    add_custom_command(
      TARGET SUNDIALS-install
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/lib/steps/"
      COMMAND
        ${CMAKE_COMMAND} -E copy
        "${tpdir}/SUNDIALS-install/lib/libsundials_${libname}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "${tpdir}/SUNDIALS-install/lib/libsundials_${libname}.6${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "${CMAKE_BINARY_DIR}/lib/steps/")
  else()
    add_custom_command(
      TARGET SUNDIALS-install
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/lib/steps/"
      COMMAND
        ${CMAKE_COMMAND} -E copy
        "${tpdir}/SUNDIALS-install/lib/libsundials_${libname}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "${tpdir}/SUNDIALS-install/lib/libsundials_${libname}${CMAKE_SHARED_LIBRARY_SUFFIX}.6"
        "${CMAKE_BINARY_DIR}/lib/steps/")
  endif()
endforeach()

# add target
include(ImportModernLib)
add_library_target(
  NAME SUNDIALS
  INCLUDE_DIRECTORIES "${tpdir}/SUNDIALS-install/include"
  LIBRARIES "${CMAKE_BINARY_DIR}/lib/steps/libsundials_nvecserial${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${CMAKE_BINARY_DIR}/lib/steps/libsundials_cvode${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${CMAKE_BINARY_DIR}/lib/steps/libsundials_ida${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${CMAKE_BINARY_DIR}/lib/steps/libsundials_kinsol${CMAKE_SHARED_LIBRARY_SUFFIX}")

# ensure that SUNDIALS is compiled and installed before the target is linked
add_dependencies(SUNDIALS::SUNDIALS SUNDIALS-install)

if(CMAKE_VERSION VERSION_LESS 3.12.0)
  add_definitions(-DSTEPS_SUNDIALS_VERSION_MAJOR=6)
else()
  add_compile_definitions(STEPS_SUNDIALS_VERSION_MAJOR=6)
endif()
set(STEPS_SUNDIALS_VERSION_MAJOR 6)
