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

// Rationale:
// We need a function to initialize stuffs such as logs
// before other STEPS object being created.
// The main function may also be used to fetch arguments
// passed from Python command.

// STEPS headers.
#include "util/common.h"

#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace mpi {

void mpiInit();
int getRank(MPI_Comm comm = MPI_COMM_WORLD);
int getNHosts(MPI_Comm comm = MPI_COMM_WORLD);

} // namespace mpi
} // namespace steps
