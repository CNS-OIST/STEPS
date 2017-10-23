# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2014 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

""" Unit tests for directional dcst."""

from __future__ import print_function
import unittest2
import random
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solv
from steps.utilities import meshio
import time

class TetDirectionalDcstTestCase(unittest2.TestCase):
    """ Tests for tetrahedral directional dcst. """
    def setUp(self):
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)
        volsys = smodel.Volsys('vsys',self.model)
        D_a = smodel.Diff('D_a', volsys, A)
        self.DCST = 0.2e-9
        D_a.setDcst(self.DCST)

        self.mesh = meshio.importAbaqus2("directional_dcst_test/mesh_tet.inp", "directional_dcst_test/mesh_tri.inp", 1e-6, "directional_dcst_test/mesh_conf")[0]
        comp = sgeom.TmComp("comp", self.mesh, range(self.mesh.ntets))
        comp.addVolsys("vsys")
    
        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
    
        self.solver = solv.Tetexact(self.model, self.mesh, self.rng)
    
        boundary_tris = self.mesh.getROIData("boundary")
        boundary_tets1 = self.mesh.getROIData("boundary_tets_1")
        boundary_tets2 = self.mesh.getROIData("boundary_tets_2")
        
        self.pairing = {}
        for tri in boundary_tris:
            neigh_tets = self.mesh.getTriTetNeighb(tri)
            if neigh_tets[0] in boundary_tets1:
                self.pairing[tri] = (neigh_tets[0], neigh_tets[1])
            else:
                self.pairing[tri] = (neigh_tets[1], neigh_tets[0])
    
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None
        self.solver = None
    
    def testTetZeroDirectionDiffRate(self):
        print("Testing tet zero directional diffusion rate from v1 to v2...")
        for tri in self.pairing.keys():
            self.solver.setTetDiffD(self.pairing[tri][0], "D_a", 0, self.pairing[tri][1])
            self.assertEqual(self.solver.getTetDiffD(self.pairing[tri][0], "D_a", self.pairing[tri][1]), 0)
            self.solver.setTetCount(self.pairing[tri][0], "A", 10)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))

    def testTetNonZeroDirectionDiffRate(self):
        print("Testing tet non-zero directional diffusion rate from v1 to v2...")
        for tri in self.pairing.keys():
            self.solver.setTetDiffD(self.pairing[tri][0], "D_a", self.DCST / 10, self.pairing[tri][1])
            self.assertEqual(self.solver.getTetDiffD(self.pairing[tri][0], "D_a", self.pairing[tri][1]), self.DCST / 10)
            self.solver.setTetCount(self.pairing[tri][0], "A", 10)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))

    def testTetDirectionalDcstCheckpoint(self):
        print("Testing checkpoint for tet directional dcst...")
        for tri in self.pairing.keys():
            self.solver.setTetDiffD(self.pairing[tri][0], "D_a", self.DCST / 10, self.pairing[tri][1])
            self.assertEqual(self.solver.getTetDiffD(self.pairing[tri][0], "D_a", self.pairing[tri][1]), self.DCST / 10)
            self.solver.setTetCount(self.pairing[tri][0], "A", 10)
        v1_count = self.solver.getROICount("v1_tets", "A")
        v2_count = self.solver.getROICount("v2_tets", "A")
        print("V1 Count: ", v1_count)
        print("V2 Count: ", v2_count)
        self.solver.checkpoint("tet_dir_dcst.cp")
        self.solver.run(1)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))
        self.solver.restore("tet_dir_dcst.cp")
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))
        self.assertEqual(self.solver.getROICount("v1_tets", "A"), v1_count)
        self.assertEqual(self.solver.getROICount("v2_tets", "A"), v2_count)

    def testTetNonDirectionDiffRate(self):
        print("Testing tet non directional diffusion rate...")
        for tri in self.pairing.keys():
            self.solver.setTetDiffD(self.pairing[tri][0], "D_a", self.DCST)
            self.assertEqual(self.solver.getTetDiffD(self.pairing[tri][0], "D_a"), self.DCST)
            self.solver.setTetCount(self.pairing[tri][0], "A", 10)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getROICount("v1_tets", "A"))
        print("V2 Count: ", self.solver.getROICount("v2_tets", "A"))

class TriDirectionalDcstTestCase(unittest2.TestCase):
    """ Tests for triangular directional dcst. """
    def setUp(self):
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)
        surfsys = smodel.Surfsys('ssys',self.model)
        D_a = smodel.Diff('D_a', surfsys, A)
        self.DCST = 0.2e-9
        D_a.setDcst(self.DCST)
        
        self.mesh = meshio.importAbaqus2("directional_dcst_test/mesh_tet.inp", "directional_dcst_test/mesh_tri.inp", 1e-6, "directional_dcst_test/mesh_conf")[0]

        boundary_tris = self.mesh.getROIData("boundary")
        v1_tets = self.mesh.getROIData("v1_tets")

        comp1 = sgeom.TmComp("comp1", self.mesh, v1_tets)
        
        patch1 = sgeom.TmPatch("patch", self.mesh, boundary_tris, comp1)
        patch1.addSurfsys("ssys")
        
        self.neigh_tris = self.mesh.getROIData("neigh_tri")
        self.focus_tri = self.mesh.getROIData("focus_tri")[0]
        
        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
        self.solver = solv.Tetexact(self.model, self.mesh, self.rng)
    
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None
        self.solver = None
    
    def testTriZeroDirectionDiffRate(self):
        print("Testing tri zero directional diffusion rate...")
        for tri in self.neigh_tris:
            self.solver.setTriSDiffD(self.focus_tri, "D_a", 0, tri)
            self.assertEqual(self.solver.getTriSDiffD(self.focus_tri, "D_a", tri), 0)
        self.solver.setTriCount(self.focus_tri, "A", 10)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
        self.solver.run(1)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
    
    def testTriNonZeroDirectionDiffRate(self):
        print("Testing tri non-zero directional diffusion rate...")
        for tri in self.neigh_tris:
            self.solver.setTriSDiffD(self.focus_tri, "D_a", self.DCST / 10, tri)
            self.assertEqual(self.solver.getTriSDiffD(self.focus_tri, "D_a", tri), self.DCST / 10)
        self.solver.setTriCount(self.focus_tri, "A", 10)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
        self.solver.run(1)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
    
    def testTriDirectionalDcstCheckpoint(self):
        print("Testing checkpoint for tri directional dcst...")
        for tri in self.neigh_tris:
            self.solver.setTriSDiffD(self.focus_tri, "D_a", self.DCST / 10, tri)
            self.assertEqual(self.solver.getTriSDiffD(self.focus_tri, "D_a", tri), self.DCST / 10)
        self.solver.setTriCount(self.focus_tri, "A", 10)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
        self.solver.checkpoint("tri_dir_dcst.cp")
        self.solver.run(1)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
        self.solver.restore("tri_dir_dcst.cp")
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
        self.assertEqual(self.solver.getPatchCount("patch", "A"), 10)
        self.assertEqual(self.solver.getTriCount(self.focus_tri, "A"), 10)

    def testTriNonDirectionDiffRate(self):
        print("Testing tri non-directional dcst...")
        self.solver.setTriSDiffD(self.focus_tri, "D_a", self.DCST)
        self.assertEqual(self.solver.getTriSDiffD(self.focus_tri, "D_a"), self.DCST)
        self.solver.setTriCount(self.focus_tri, "A", 10)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))
        self.solver.run(1)
        print("Patch Count: ", self.solver.getPatchCount("patch", "A"))
        print("tri Count: ", self.solver.getTriCount(self.focus_tri, "A"))

class DiffBndDirectionalDcstTestCase(unittest2.TestCase):
    """ Tests for directional dcst in diffusion boundary. """
    def setUp(self):
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)
        volsys = smodel.Volsys('vsys',self.model)
        D_a = smodel.Diff('D_a', volsys, A)
        self.DCST = 0.2e-9
        D_a.setDcst(self.DCST)
        
        self.mesh = meshio.importAbaqus2("directional_dcst_test/mesh_tet.inp", "directional_dcst_test/mesh_tri.inp", 1e-6, "directional_dcst_test/mesh_conf")[0]

        boundary_tris = self.mesh.getROIData("boundary")
        v1_tets = self.mesh.getROIData("v1_tets")
        v2_tets = self.mesh.getROIData("v2_tets")

        comp1 = sgeom.TmComp("comp1", self.mesh, v1_tets)
        comp2 = sgeom.TmComp("comp2", self.mesh, v2_tets)
        
        comp1.addVolsys("vsys")
        comp2.addVolsys("vsys")
        
        db = sgeom.DiffBoundary("boundary", self.mesh, boundary_tris)
        
        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
        self.solver = solv.Tetexact(self.model, self.mesh, self.rng)
    
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None
        self.solver = None
    
    def testDBZeroDirectionDiffRate(self):
        print("Testing diff bnd zero directional diffusion rate...")
        self.solver.setCompCount("comp1", "A", 100)
        self.solver.setDiffBoundaryDcst("boundary", "A", 0.0)
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))
    
    def testBDNonZeroDirectionDiffRate(self):
        print("Testing diff bnd non-zero directional diffusion rate")
        self.solver.setCompCount("comp1", "A", 100)
        self.solver.setDiffBoundaryDcst("boundary", "A", self.DCST / 10, "comp2")
        self.solver.setDiffBoundaryDcst("boundary", "A", 0.0, "comp1")
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))
    
    def testDBSingleDirectionDiffRate(self):
        print("Testing diff bnd nonzero single direction diffusion rate...")
        self.solver.setCompCount("comp1", "A", 100)
        self.solver.setDiffBoundaryDcst("boundary", "A", self.DCST, "comp2" )
        self.solver.setDiffBoundaryDcst("boundary", "A" , 0.0, "comp1")
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))

    def testNonzeroBidirectionDiffRate(self):
        print("Testing diff bnd nonzero bidirection diffusion rate...")
        self.solver.setCompCount("comp1", "A", 100)
        self.solver.setCompCount("comp2", "A", 10)
        self.solver.setDiffBoundaryDcst("boundary", "A", self.DCST / 10)
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))
        self.solver.run(1)
        print("V1 Count: ", self.solver.getCompCount("comp1", "A"))
        print("V2 Count: ", self.solver.getCompCount("comp2", "A"))

def suite():
    all_tests = []
    all_tests.append(unittest2.makeSuite(TetDirectionalDcstTestCase, "test"))
    all_tests.append(unittest2.makeSuite(TriDirectionalDcstTestCase, "test"))
    all_tests.append(unittest2.makeSuite(DiffBndDirectionalDcstTestCase, "test"))
    return unittest2.TestSuite(all_tests)

if __name__ == "__main__":
    unittest2.TextTestRunner(verbosity=2).run(suite())
