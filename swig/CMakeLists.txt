####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


include_directories(. ../src)
include_directories(${PYTHON_INCLUDE_DIRS})

include(UseSWIG)

# UseSWIG does not do dependency analysis. Great work, CMake developers.

set(swig_included_sources
    "error.i"    "geom.i"    "model.i"    "rng.i"    "solver.i")

if(NUMPY_FOUND)
    set(swig_sources "steps_swig_numpy.i")
    include_directories(${NUMPY_INCLUDE_DIRS})
    set(swig_module "steps_swig_numpy")

    list(APPEND swig_included_sources "numpy.i" "unchecked_stl_seq.i")
else()
    set(swig_sources "steps_swig.i")
    set(swig_module "steps_swig")
endif()


set_source_files_properties(${swig_sources} PROPERTIES CPLUSPLUS ON)

# if it is possible, give WITH_MPI def to swig if MPI_FOUND so that it can choose to compile MPI modules or not.

if(MPI_FOUND)
    set_source_files_properties(${swig_sources} PROPERTIES SWIG_FLAGS "-DWITH_MPI")
    list(APPEND swig_included_sources "mpi_solver.i")
endif()

set(SWIG_MODULE_${swig_module}_EXTRA_DEPS ${swig_included_sources})

swig_add_module(${swig_module} python ${swig_sources} $<TARGET_OBJECTS:libsteps-obj>)

set(_swig_link_libs ${libsteps_link_libraries})
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Linker flag to Skip unresolved symbols
    list(APPEND _swig_link_libs "-Wl,-undefined,dynamic_lookup")
endif()
swig_link_libraries(${swig_module} ${_swig_link_libs})


# set steps_swig_target and steps_swig_pymodules for use by py/CMakeLists.txt

set(steps_swig_target ${SWIG_MODULE_${swig_module}_REAL_NAME} PARENT_SCOPE)

set(pymodules)
foreach(ifile IN LISTS swig_sources)
    string(REGEX REPLACE  "\\.i$" ".py" pyfile ${ifile})
    list(APPEND pymodules "${CMAKE_CURRENT_BINARY_DIR}/${pyfile}")
endforeach()

set(steps_swig_pymodules ${pymodules} PARENT_SCOPE)
