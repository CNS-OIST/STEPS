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
import steps.mpi
import steps.mpi.solver as solv
from steps.utilities import meshio
import steps.utilities.geom_decompose as gd

class ParallelOpSplitTestCase(unittest.TestCase):
    """ Test cases for parallel OpSplit solver. """
    def setUp(self):
        self.model = smodel.Model()
        A = smodel.Spec("A", self.model)
        B = smodel.Spec("B", self.model)
        C = smodel.Spec("C", self.model)
        D = smodel.Spec("D", self.model)
        E = smodel.Spec("E", self.model)

        self.vsys1 = smodel.Volsys('vsys1', self.model)
        self.vsys2 = smodel.Volsys('vsys2', self.model)
        self.ssys1 = smodel.Surfsys('ssys1', self.model)

        self.reac1 = smodel.Reac('reac1', self.vsys1, lhs = [A], rhs = [A],  kcst = 1e5)
        self.reac2 = smodel.Reac('reac2', self.vsys2, lhs = [E], rhs = [E],  kcst = 1e4)

        self.sreac = smodel.SReac('sreac', self.ssys1, slhs = [B], srhs = [B],  kcst = 1e3)
    
        if __name__ == "__main__":
            self.mesh = meshio.loadMesh('meshes/cyl_len10_diam1')[0]
        else:
            self.mesh = meshio.loadMesh('getROIArea_bugfix_test/meshes/cyl_len10_diam1')[0]

        ntets = self.mesh.countTets()
        comp1Tets, comp2Tets = [], []
        comp1Tris, comp2Tris = set(), set()
        for i in range(ntets):
            if self.mesh.getTetBarycenter(i)[0] > 0:
                comp1Tets.append(i)
                comp1Tris |= set(self.mesh.getTetTriNeighb(i))
            else:
                comp2Tets.append(i)
                comp2Tris |= set(self.mesh.getTetTriNeighb(i))
        patch1Tris = list(comp1Tris & comp2Tris)

        self.comp1 = sgeom.TmComp('comp1', self.mesh, comp1Tets)
        self.comp2 = sgeom.TmComp('comp2', self.mesh, comp2Tets)
        self.comp1.addVolsys('vsys1')
        self.comp2.addVolsys('vsys2')

        self.patch1 = sgeom.TmPatch('patch1', self.mesh, patch1Tris, self.comp1, self.comp2)
        self.patch1.addSurfsys('ssys1')

        self.ROI1 = self.mesh.addROI('ROI1', sgeom.ELEM_TET, comp1Tets)
        self.ROI2 = self.mesh.addROI('ROI2', sgeom.ELEM_TET, comp2Tets)
        self.ROI3 = self.mesh.addROI('ROI3', sgeom.ELEM_TRI, patch1Tris)
    
        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        tet_hosts = gd.linearPartition(self.mesh, [steps.mpi.nhosts, 1, 1])
        tri_hosts = gd.partitionTris(self.mesh, tet_hosts, patch1Tris)
        self.solver = solv.TetOpSplit(self.model, self.mesh, self.rng, solv.EF_NONE, tet_hosts, tri_hosts)

        
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None
        self.solver = None

    def testsetROICount(self):
        """
        test for bugfix https://github.com/CNS-OIST/HBP_STEPS/issues/259
        """
        self.solver.setROICount("ROI1", "A", 100)
        get_count = self.solver.getROICount("ROI1", 'A')
        self.assertEqual(get_count, 100)
        h_mu = self.solver.getCompReacA("comp1", "reac1")
        self.assertNotEqual(h_mu, 0.0)

        self.solver.setROICount("ROI2", "E", 200)
        get_count = self.solver.getROICount("ROI2", 'E')
        self.assertEqual(get_count, 200)
        h_mu = self.solver.getCompReacA("comp2", "reac2")
        self.assertNotEqual(h_mu, 0.0)

        self.solver.setROICount("ROI3", "B", 300)
        get_count = self.solver.getROICount("ROI3", 'B')
        self.assertEqual(get_count, 300)
        h_mu = self.solver.getPatchSReacA("patch1", "sreac")
        self.assertNotEqual(h_mu, 0.0)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ParallelOpSplitTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
