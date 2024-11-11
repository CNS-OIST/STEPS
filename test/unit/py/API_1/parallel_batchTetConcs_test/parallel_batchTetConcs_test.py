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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.mpi as smpi
import steps.mpi.solver as spsolver
import steps.utilities.geom_decompose as sgd
import numpy as np

from steps.API_1.geom import INDEX_DTYPE

class batchTetConcs(unittest.TestCase):
    """
    Test batch set/getBatchTetConcs both list and np.array version
    """
    def setUp(self):
        # model setup
        mdl = smodel.Model()
        spc = smodel.Spec('A', mdl)
        vsys = smodel.Volsys('A', mdl)
        diff = smodel.Diff('diff_A', vsys, spc)
        diff.setDcst(0.0)

        # mesh
        vertCoos = [0.0, 0.0, 0.0, \
                    1.0e-6, 0.0, 0.0, \
                    0.0, 1.0e-6, 0.0, \
                    0.0, 0.0, 1.0e-6, \
                    1.0e-6, 1.0e-6, 1.0e-6 ]
        vertIds = [0, 1, 2, 3, \
                   1, 2, 3, 4  ]

        # geom setup
        msh = sgeom.Tetmesh(vertCoos, vertIds)
        ntets = msh.countTets()
        comp = sgeom.TmComp('comp', msh, range(ntets))
        comp.addVolsys('A')

        # init sim
        rng = srng.create('mt19937', 512)
        rng.initialize(2903)

        tet_hosts = sgd.linearPartition(msh, [1, 1, smpi.nhosts])
        self.sim = spsolver.TetOpSplit(mdl, msh, rng, spsolver.EF_NONE, tet_hosts)

    def testSetBatchTetConcs(self):
        self.sim.reset()
        tetIds = [0, 1]
        tetConcs = [1.0, 2.0]

        # set
        self.sim.setBatchTetConcs(tetIds, 'A', tetConcs)

        # check
        for i in range(0, len(tetIds)):
            self.assertAlmostEqual(self.sim.getTetSpecConc(tetIds[i], 'A'), tetConcs[i])

    @unittest.expectedFailure
    def testSetBatchTetConcsSizeFail(self):
        self.sim.reset()
        tetIds = [0, 1, 2]
        tetConcs = [1.0, 2.0]

        # set
        self.sim.setBatchTetConcs(tetIds, 'A', tetConcs)

    @unittest.expectedFailure
    def testSetBatchTetConcsOutOfRangeFail(self):
       self.sim.reset()
       tetIds = [0, 10]
       tetConcs = [1.0, 2.0]

       # set
       self.sim.setBatchTetConcs(tetIds, 'A', tetConcs)

    def testGetBatchTetConcs(self):
        self.sim.reset()
        tetIds = [0, 1]
        tetConcs = [1.0, 2.0]

        # set
        for i in range(0, len(tetIds)):
            self.sim.setTetSpecConc(tetIds[i], 'A', tetConcs[i])

        # get
        tetConcsBatch = self.sim.getBatchTetConcs(tetIds, 'A')

    @unittest.expectedFailure
    def testGetBatchTetConcsOutOfRangeFail(self):
       self.sim.reset()
       tetIds = [0, 10]

       # get
       tetConcs = self.sim.getBatchTetConcs(tetIds, 'A')

    def testSetBatchTetSpecConcsNP(self):
        self.sim.reset()
        tetIds = np.array([0, 1], dtype=INDEX_DTYPE)
        tetConcs = np.array([1.0, 2.0], dtype=float)

        # set
        self.sim.setBatchTetSpecConcsNP(tetIds, 'A', tetConcs)

        # check
        for i in range(0, len(tetIds)):
            self.assertAlmostEqual(self.sim.getTetSpecConc(tetIds[i], 'A'), tetConcs[i])

    def testGetBatchTetConcsNP(self):
        self.sim.reset()
        tetIds = np.array([0, 1], dtype=INDEX_DTYPE)
        tetConcs = np.array([1.0, 2.0], dtype=float)

        # set
        for i in range(0, len(tetIds)):
            self.sim.setTetSpecConc(tetIds[i], 'A', tetConcs[i])

        # get
        tetConcsBatch = np.zeros((len(tetConcs),), dtype=float)
        self.sim.getBatchTetConcsNP(tetIds, 'A', tetConcsBatch)

        # check
        for i in range(0, len(tetIds)):
            self.assertAlmostEqual(tetConcsBatch[i], tetConcs[i])

    @unittest.expectedFailure
    def testGetBatchTetConcsNPSizeFail(self):
        self.sim.reset()
        tetIds = np.array([0, 1], dtype=INDEX_DTYPE)
        tetConcs = np.array([1.0, 2.0], dtype=float)

        # set
        for i in range(0, len(tetIds)):
            self.sim.setTetSpecConc(tetIds[i], 'A', tetConcs[i])

        # get
        tetConcsBatch = np.zeros((len(tetConcs)+1,), dtype=float)
        self.sim.getBatchTetConcsNP(tetIds, 'A', tetConcsBatch)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(batchTetConcs))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
