/*
 #################################################################################
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

 */

#include "finish.hpp"

#ifdef OMEGA_H_USE_GMSH
#include <gmsh.h>
#endif // OMEGA_H_USE_GMSH

#include "util/profile/profiler_interface.h"

void steps::finish() {
#ifdef OMEGA_H_USE_GMSH
    ::gmsh::finalize();
#endif // OMEGA_H_USE_GMSH

    steps::Instrumentor::finalize_profile();
}
