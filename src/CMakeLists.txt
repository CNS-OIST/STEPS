# #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

if(STEPS_ENABLE_ERROR_ON_WARNING)
  add_options(ALL_LANGUAGES ALL_BUILDS -Werror)
endif()

add_definitions(-DENABLE_ASSERTLOG=1)

set(lib_soversion "${VERSION_MAJOR}.${VERSION_MINOR}")
set(lib_version "${lib_soversion}.${VERSION_PATCH}")

add_subdirectory(steps/math)
add_subdirectory(steps/geom)
add_subdirectory(steps/model)
add_subdirectory(steps/rng)
add_subdirectory(steps/solver)
add_subdirectory(steps/tetexact)
add_subdirectory(steps/tetode)
add_subdirectory(steps/util)
add_subdirectory(steps/wmdirect)
add_subdirectory(steps/wmrk4)
add_subdirectory(steps/wmrssa)

# enable below to turn on MPI profiling
if(USE_MPI)
  add_subdirectory(steps/mpi)
  if(USE_MPI_PROFILING)
    add_definitions(-DMPI_PROFILING=1)
  endif()
endif()

# ==============================================================================
# ==============================

add_library(libsteps_static STATIC $<TARGET_OBJECTS:stepswmrssa>)
set_target_properties(libsteps_static PROPERTIES OUTPUT_NAME steps)
target_link_libraries(libsteps_static PUBLIC
  stepsutil
  stepsmath
  stepsgeom
  stepsrng
  stepsmodel
  stepssolver
  stepstetexact
  stepstetode
  stepswmdirect
  stepswmrk4
  ${libsteps_link_libraries}
)
if(USE_MPI)
  target_link_libraries(libsteps_static PUBLIC stepstetopsplit)
  if(STEPS_USE_DIST_MESH)
    target_link_libraries(libsteps_static PUBLIC stepsdist stepsdisttest)
  endif()
endif()

add_library(libsteps SHARED $<TARGET_OBJECTS:stepswmrssa>)
target_link_libraries(libsteps PUBLIC
  stepsutil
  stepsmath
  stepsgeom
  stepsrng
  stepsmodel
  stepssolver
  stepstetexact
  stepstetode
  stepswmdirect
  stepswmrk4
  ${libsteps_link_libraries}
)
if(USE_MPI)
  target_link_libraries(libsteps PUBLIC stepstetopsplit)
  if(STEPS_USE_DIST_MESH)
    target_link_libraries(libsteps PUBLIC stepsdist stepsdisttest)
  endif()
endif()

set_target_properties(libsteps
                      PROPERTIES VERSION
                                 "${lib_version}"
                                 SOVERSION
                                 "${lib_soversion}")
set_target_properties(libsteps PROPERTIES OUTPUT_NAME steps)

if(STEPS_USE_DIST_MESH)
    add_executable(tetopsplit_dist steps/mpi/dist/tetopsplit/client/main.cpp)
    target_link_libraries(tetopsplit_dist PUBLIC libsteps stepsdisttest)
endif()
# ==========================================================================
# Installation of the headers and shared lib
#
option(INSTALL_DEV_LIBS
       "Install steps lib and headers which are required for development only"
       OFF)

if(INSTALL_DEV_LIBS)
  install(TARGETS libsteps
          RUNTIME DESTINATION bin COMPONENT lib
          ARCHIVE DESTINATION lib COMPONENT lib
          LIBRARY DESTINATION lib COMPONENT lib)
endif()
