# -*- coding: utf-8 -*-

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
from __future__ import print_function

try:
    from .. import cysteps_mpi as stepslib
except Exception as e:
    raise Exception("[ERROR] Could not load cysteps_mpi.so. Please verify it exists and system steps was built with MPI support.\n Details: " + str(e))

stepslib.mpiInit()
import atexit
atexit.register(stepslib.mpiFinish)

rank = stepslib.getRank()
nhosts = stepslib.getNHosts()

import steps
if not steps._quiet and rank == 0:
    print("-----------------------------------------------------------------")
    print("STEPS is running in parallel mode with ", nhosts, " processes")
    print("-----------------------------------------------------------------")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
