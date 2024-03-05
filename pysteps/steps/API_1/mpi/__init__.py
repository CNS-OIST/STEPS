# -*- coding: utf-8 -*-

####################################################################################
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
from __future__ import print_function

import sys

from steps import stepslib

_mpiSuffixes = ['mpi', 'dist']
if not any(stepslib.__name__.endswith(suff) for suff in _mpiSuffixes):
    raise ImportError(
        f'[ERROR] Could not load cysteps_mpi.so. Please verify it exists and system steps was '
        f'built with MPI support.'
    )

# Force stderr flush when an exception is raised.
# Without this, under some conditions, if one process raises a python exception,
# the corresponding message is not always printed out as it should be.
def customHook(tpe, val, bt):
    sys.__excepthook__(tpe, val, bt)
    sys.stderr.flush()
    stepslib.mpiAbort()

sys.excepthook = customHook

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
