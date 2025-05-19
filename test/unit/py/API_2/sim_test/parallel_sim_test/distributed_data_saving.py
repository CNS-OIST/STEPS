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

""" Unit tests for distributed data saving."""

import unittest
import sys

from steps import interface

from steps.geom import *
from steps.rng import *
from steps.sim import *

from .. import test_dataSaving as tds
from .. import base_model as bm


class DistTetOpSplitSimDataSaving(bm.DistTetTestModelFramework, tds.TetSimDataSaving):
    """Test data access, setting, and saving with tetmeshes."""

    def setUp(self):
        bm.DistTetTestModelFramework.setUp(self)
        tds.TetSimDataSaving.setUp(self, False)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('DistTetOpSplit', nmdl, ngeom, nrng)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        return None

    @unittest.skip('Not needed here')
    def testXDMFWithoutMPI(self):
        pass


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistTetOpSplitSimDataSaving))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
