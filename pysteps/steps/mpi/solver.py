# -*- coding: utf-8 -*-

####################################################################################
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
###

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# This file is the user-interface file for all solver objects.
# All objects are directly derived from the corresponding swig objects.
# All objects are owned by Python.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
"""
Implementation of parallel simulation solvers.

Each solver is a partial or full implementation of the STEPS solver API.

The steps.mpi.solver.TetOpSplit class.

"""
from steps import stepslib
from steps.solver import _Base_Solver

# Constants aliases (yep, must be hand coded)
EF_NONE = stepslib._py_API.EF_NONE
EF_DEFAULT = stepslib._py_API.EF_DEFAULT
EF_DV_BDSYS = stepslib._py_API.EF_DV_BDSYS
EF_DV_SLUSYS = stepslib._py_API.EF_DV_SLUSYS
EF_DV_PETSC  = stepslib._py_API.EF_DV_PETSC

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Tetrahedral Direct SSA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #        
class TetOpSplit(stepslib._py_TetOpSplitP, _Base_Solver) :
    """
    Construction::
    
        sim = steps.solver.TetOpSplit(model, geom, rng, tet_hosts=[], tri_hosts={}, wm_hosts=[], calcMembPot=0)
    
    Create a spatial stochastic solver based on operator splitting, that is that reaction events are partitioned and diffusion is approximated. 
    If voltage is to be simulated, argument calcMembPot specifies the solver e.g. calcMembPot=steps.solver.EF_DV_PETSC will utilise the PETSc library. calcMembPot=0 means voltage will not be simulated. 
    
    Arguments:
    steps.model.Model model
    steps.geom.Geom geom
    steps.rng.RNG rng
    list<int> tet_hosts (default=[])
    dict<int, int> tri_hosts (default={})
    list<int> wm_hosts (default=[])
    int calcMemPot (default=0)
    
    """
    def run(self, end_time, cp_interval=0.0, prefix=""):
        """
        Run the simulation until <end_time>, 
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'tetopsplitP' )
        
    def advance(self, advance_time, cp_interval=0.0, prefix=""):
        """
        Advance the simulation for <advance_time>, 
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        end_time = self.getTime() + advance_time
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'tetopsplitP')

    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        return self._getIndexMapping()
