/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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

 */


#include <mpi.h>

#include "steps/mpi/mpi_finish.hpp"
#include "steps/mpi/mpi_common.hpp"

#include "easylogging++.h"

#ifdef USE_PETSC
#include "petscsys.h"
#endif

void steps::mpi::mpiFinish(void) {
#ifdef USE_PETSC
    PetscBool PETSc_used = PETSC_FALSE;
    PetscInitialized(&PETSc_used);

    if (!PETSc_used) {
    int status = MPI_Finalize();
    }
#endif
}

////////////////////////////////////////////////////////////////////////////////
