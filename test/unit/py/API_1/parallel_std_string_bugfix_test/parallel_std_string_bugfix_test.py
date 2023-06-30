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

# Test if diffusion selector works as expected
# Used to detect bug in https://github.com/CNS-OIST/HBP_STEPS/issues/157

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.mpi
import steps.mpi.solver as solv
from steps.utilities import meshio
import steps.utilities.geom_decompose as gd

import numpy as np

class ParallelStdStringBugfixCase(unittest.TestCase):
    """ 
    Test if cython bindings for TetOpSplit work as expected.
    The folowing methods did not use to_std_string for converting the python
    strings to std strings:
        checkpoint
        restore
        saveMembOpt
    o   getBatchTetCounts
    o   getBatchTriCounts
    o   setBatchTetConcs
    o   getBatchTetConcs
    o   getBatchTetCountsNP
    o   getBatchTriCountsNP
    o   sumBatchTetCountsNP
    o   sumBatchTriCountsNP
        sumBatchTriGHKIsNP
        sumBatchTriOhmicIsNP
    The tests are covering all methods (marked by an o) that do not require PETSC.
    checkpoint and restore are not implemented in TetOpSplit so they cannot be 
    tested either.
    """
    def setUp(self):
        DCST = 0.08e-10
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)

        self.vsys = smodel.Volsys('vsys', self.model)
        self.ssys = smodel.Surfsys('ssys', self.model)
        self.diff = smodel.Diff("diff", self.vsys, A, DCST)
        self.sdiff = smodel.Diff("diff", self.ssys, A, DCST)

        if __name__ == "__main__":
            self.mesh = meshio.importAbaqus('meshes/test_mesh.inp', 1e-7)[0]
        else:
            self.mesh = meshio.importAbaqus('parallel_std_string_bugfix_test/meshes/test_mesh.inp', 1e-7)[0]

        self.tmcomp = sgeom.TmComp('comp', self.mesh, range(self.mesh.ntets))
        self.tmcomp.addVolsys('vsys')
        self.surf_tris = self.mesh.getSurfTris()
        self.tmpatch = sgeom.TmPatch('patch', self.mesh, self.surf_tris, icomp = self.tmcomp)
        self.tmpatch.addSurfsys('ssys')

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
        tet_hosts = gd.binTetsByAxis(self.mesh, steps.mpi.nhosts)
        tri_hosts = gd.partitionTris(self.mesh, tet_hosts, self.surf_tris)
        self.solver = solv.TetOpSplit(self.model, self.mesh, self.rng, solv.EF_NONE, tet_hosts, tri_hosts)

        self.solver.reset()

        
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None

    def testGetBatchTetCounts(self):
        tetCounts = [i % 15 for i in range(self.mesh.ntets)]
        for i, v in enumerate(tetCounts):
            self.solver.setTetSpecCount(i, 'A', v)

        tetCounts2 = self.solver.getBatchTetSpecCounts(range(self.mesh.ntets), 'A')
        self.assertListEqual(tetCounts, tetCounts2)

    def testGetBatchTriCounts(self):
        triCounts = [i % 20 for i in self.surf_tris]
        for i, v in zip(self.surf_tris, triCounts):
            self.solver.setTriSpecCount(i, 'A', v)
        triCounts2 = self.solver.getBatchTriSpecCounts(self.surf_tris, 'A')
        self.assertListEqual(triCounts, triCounts2)

    def testSetGetBatchTetConcs(self):
        Nav = 6.02214076e23
        tetConcs = [(i % 15) / Nav / self.mesh.getTetVol(i) for i in range(self.mesh.ntets)]
        self.solver.setBatchTetConcs(range(self.mesh.ntets), 'A', tetConcs)

        tetConcs2 = self.solver.getBatchTetConcs(range(self.mesh.ntets), 'A')

        for c1, c2 in zip(tetConcs, tetConcs2):
            self.assertLessEqual(abs(c1 - c2), 0.1 / 100 * c1)

    def testGetBatchTetCountsNP(self):
        tetCounts = [i % 15 for i in range(self.mesh.ntets)]
        for i, v in enumerate(tetCounts):
            self.solver.setTetSpecCount(i, 'A', v)

        tetInds = np.array(range(self.mesh.ntets), dtype=sgeom.INDEX_DTYPE)
        tetCounts2 = np.zeros(self.mesh.ntets)
        self.solver.getBatchTetSpecCountsNP(tetInds, 'A', tetCounts2)

        for c1, c2 in zip(tetCounts, tetCounts2):
            self.assertEqual(c1, c2)

    def testGetBatchTriCountsNP(self):
        triCounts = [i % 20 for i in self.surf_tris]
        for i, v in zip(self.surf_tris, triCounts):
            self.solver.setTriSpecCount(i, 'A', v)

        triInds = np.array(self.surf_tris, dtype=sgeom.INDEX_DTYPE)
        triCounts2 = np.zeros(len(self.surf_tris))
        self.solver.getBatchTriSpecCountsNP(triInds, 'A', triCounts2)

        for c1, c2 in zip(triCounts, triCounts2):
            self.assertEqual(c1, c2)

    def testSumMatchTetCountsNP(self):
        tetCounts = [i % 15 for i in range(self.mesh.ntets)]
        for i, v in enumerate(tetCounts):
            self.solver.setTetSpecCount(i, 'A', v)

        tetInds = np.array(range(self.mesh.ntets), dtype=sgeom.INDEX_DTYPE)
        tetSum = self.solver.sumBatchTetCountsNP(tetInds, 'A')

        self.assertEqual(tetSum, sum(tetCounts))

    def testSumMatchTriCountsNP(self):
        triCounts = [i % 20 for i in self.surf_tris]
        for i, v in zip(self.surf_tris, triCounts):
            self.solver.setTriSpecCount(i, 'A', v)

        triInds = np.array(self.surf_tris, dtype=sgeom.INDEX_DTYPE)
        triSum = self.solver.sumBatchTriCountsNP(triInds, 'A')

        self.assertEqual(triSum, sum(triCounts))


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ParallelStdStringBugfixCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

