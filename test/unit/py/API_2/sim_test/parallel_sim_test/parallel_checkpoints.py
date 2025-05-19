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

""" Unit tests for TetOpSplit and TetVesicle checkpointing."""

import importlib
import unittest

from steps import interface

from steps.geom import *
from steps.rng import *
from steps.sim import *

from .. import base_model
from .. import test_checkpoints


class TetOpSplitCheckpoints(test_checkpoints.TetCheckpoints):
    """Test checkpointing for TetOpSplit"""

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        part = LinearMeshPartition(ngeom, 1, MPI.nhosts, 1)
        nsim = Simulation('TetOpSplit', nmdl, ngeom, nrng, MPI.EF_DEFAULT, part)
        nsim.EfieldDT = self.efielddt
        return nsim


class TetVesicleCheckpoints(test_checkpoints.TetCheckpoints, base_model.VesTestModelFramework):
    """Test checkpointing for TetVesicle"""

    def setUp(self):
        base_model.VesTestModelFramework.setUp(self)
        self.useEField = True
        self.useVesicle = True

        test_checkpoints.WmCheckpoints.setUp(self, False)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('TetVesicle', nmdl, ngeom, nrng, MPI.EF_DEFAULT)
        nsim.EfieldDT = self.efielddt
        self.one_time_init_API2_sim(nsim)
        return nsim

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetOpSplitCheckpoints))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetVesicleCheckpoints))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


