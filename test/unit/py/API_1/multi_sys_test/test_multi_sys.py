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

# Provide tests for more complicate volume/surface systems, as current validations
# mostly focus on single vsys/ssys, in which case the validation can not identify multi system issues.
# For example global indices being used as local indices.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solv
from steps.utilities import meshio

class MultiVolsysTestCase(unittest.TestCase):
    """ 
    Test for two volsys with completely different reactions. 
    In this case the global indices of species in the second vsys will be different from the local indices.
    """
    def setUp(self):
        KCST = 10000.0
        DCST = 0.08e-12
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)
        B = smodel.Spec('B', self.model)
        C = smodel.Spec('C', self.model)
        D = smodel.Spec('D', self.model)
        E = smodel.Spec('E', self.model)
        F = smodel.Spec('F', self.model)

        self.vsys1 = smodel.Volsys('vsys1', self.model)
        self.vsys2 = smodel.Volsys('vsys2', self.model)
    
        self.reac1 = smodel.Reac('reac1', self.vsys1, lhs = [A, B], rhs = [C],  kcst = KCST)
        self.reac2 = smodel.Reac('reac2', self.vsys2, lhs = [D, E], rhs = [F],  kcst = KCST)
        
        self.geom = sgeom.Geom()
        self.comp1 = sgeom.Comp('comp1', self.geom, 1e-18)
        self.comp1.addVolsys('vsys1')
        self.comp2 = sgeom.Comp('comp2', self.geom, 1e-18)
        self.comp2.addVolsys('vsys2')
    
        if __name__ == "__main__":
            self.mesh = meshio.importAbaqus('meshes/brick_40_4_4_1400tets.inp', 1e-6)[0]
        else:
            self.mesh = meshio.importAbaqus('multi_sys_test/meshes/brick_40_4_4_1400tets.inp', 1e-6)[0]

        comp1_tets = []
        comp2_tets = []

        for t in range(self.mesh.ntets):
            cord = self.mesh.getTetBarycenter(t)
            if cord[0] < 0.0:
                comp1_tets.append(t)
            else:
                comp2_tets.append(t)
        
        self.tmcomp1 = sgeom.TmComp('comp1', self.mesh, comp1_tets)
        self.tmcomp1.addVolsys('vsys1')
        self.tmcomp2 = sgeom.TmComp('comp2', self.mesh, comp2_tets)
        self.tmcomp2.addVolsys('vsys2')

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
    def tearDown(self):
        self.model = None
        self.geom = None
        self.mesh = None
        self.rng = None

    def _runTest(self, solver):
        solver.setCompSpecCount('comp1', 'A', 100000)
        solver.setCompSpecCount('comp1', 'B', 100000)
        solver.setCompSpecCount('comp2', 'D', 100000)
        solver.setCompSpecCount('comp2', 'E', 100000)
        solver.run(1)
        if __name__ == "__main__":
            print("")
            print("A:", solver.getCompSpecCount("comp1", "A"))
            print("B:", solver.getCompSpecCount("comp1", "B"))
            print("C:", solver.getCompSpecCount("comp1", "C"))
            print("D:", solver.getCompSpecCount("comp2", "D"))
            print("E:", solver.getCompSpecCount("comp2", "E"))
            print("F:", solver.getCompSpecCount("comp2", "F"))
        self.assertNotEqual(solver.getCompSpecCount("comp1", "C"), 0)
        self.assertNotEqual(solver.getCompSpecCount("comp2", "F"), 0)
        self.assertEqual(solver.getCompSpecCount("comp1", "A"), solver.getCompSpecCount("comp1", "B"))
        self.assertEqual(solver.getCompSpecCount("comp1", "A") + solver.getCompSpecCount("comp1", "C"), 100000)
        self.assertEqual(solver.getCompSpecCount("comp2", "D"), solver.getCompSpecCount("comp2", "E"))
        self.assertEqual(solver.getCompSpecCount("comp2", "D") + solver.getCompSpecCount("comp2", "F"), 100000)

    def testWmdirect(self):
        solver = solv.Wmdirect(self.model, self.geom, self.rng)
        self._runTest(solver)
    
    def testWmrssa(self):
        solver = solv.Wmrssa(self.model, self.geom, self.rng)
        self._runTest(solver)

    def testTetexact(self):
        solver = solv.Tetexact(self.model, self.mesh, self.rng)
        self._runTest(solver)

class MultiSurfsysTestCase(unittest.TestCase):
    """ 
    Test for two surface systems with completely different reactions. 
    In this case the global indices of species in the second ssys will be different from the local indices.
    """
    def setUp(self):
        KCST = 1e6
        DCST = 0.08e-12
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)
        B = smodel.Spec('B', self.model)
        C = smodel.Spec('C', self.model)
        D = smodel.Spec('D', self.model)
        E = smodel.Spec('E', self.model)
        F = smodel.Spec('F', self.model)

        self.ssys1 = smodel.Surfsys('ssys1', self.model)
        self.ssys2 = smodel.Surfsys('ssys2', self.model)
    
        self.sreac1 = smodel.SReac('sreac1', self.ssys1, slhs = [A, B], srhs = [C],  kcst = KCST)
        self.sreac2 = smodel.SReac('sreac2', self.ssys2, slhs = [D, E], srhs = [F],  kcst = KCST)
        
        self.geom = sgeom.Geom()
        self.comp = sgeom.Comp('comp', self.geom, 1e-18)
        self.patch1 = sgeom.Patch('patch1', self.geom, self.comp, None, 1e-12)
        self.patch1.addSurfsys('ssys1')
        self.patch2 = sgeom.Patch('patch2', self.geom, self.comp, None, 1e-12)
        self.patch2.addSurfsys('ssys2')

        if __name__ == "__main__":
            self.mesh = meshio.importAbaqus('meshes/brick_40_4_4_1400tets.inp', 1e-6)[0]
        else:
            self.mesh = meshio.importAbaqus('multi_sys_test/meshes/brick_40_4_4_1400tets.inp', 1e-6)[0]

        comp1_tets = []
        comp2_tets = []

        for t in range(self.mesh.ntets):
            cord = self.mesh.getTetBarycenter(t)
            if cord[0] < 0.0:
                comp1_tets.append(t)
            else:
                comp2_tets.append(t)
        
        self.tmcomp = sgeom.TmComp('comp', self.mesh, range(self.mesh.ntets))

        surf_tris = self.mesh.getSurfTris()

        patch_tris1 = []
        patch_tris2 = []
        for tri in surf_tris:
            tet_neighs = self.mesh.getTriTetNeighb(tri)
            for tet in tet_neighs:
                if tet in comp1_tets:
                    patch_tris1.append(tri)
                    break
                elif tet in comp2_tets:
                    patch_tris2.append(tri)
                    break
        self.tmpatch1 = sgeom.TmPatch('patch1', self.mesh, patch_tris1, self.tmcomp)
        self.tmpatch1.addSurfsys('ssys1')
        self.tmpatch2 = sgeom.TmPatch('patch2', self.mesh, patch_tris2, self.tmcomp)
        self.tmpatch2.addSurfsys('ssys2')

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
    def tearDown(self):
        self.model = None
        self.geom = None
        self.mesh = None
        self.rng = None

    def _runTest(self, solver):
        solver.setPatchSpecCount('patch1', 'A', 100000)
        solver.setPatchSpecCount('patch1', 'B', 100000)
        solver.setPatchSpecCount('patch2', 'D', 100000)
        solver.setPatchSpecCount('patch2', 'E', 100000)
        solver.run(1)
        if __name__ == "__main__":
            print("")
            print("A:", solver.getPatchSpecCount("patch1", "A"))
            print("B:", solver.getPatchSpecCount("patch1", "B"))
            print("C:", solver.getPatchSpecCount("patch1", "C"))
            print("D:", solver.getPatchSpecCount("patch2", "D"))
            print("E:", solver.getPatchSpecCount("patch2", "E"))
            print("F:", solver.getPatchSpecCount("patch2", "F"))
        self.assertNotEqual(solver.getPatchSpecCount("patch1", "C"), 0)
        self.assertNotEqual(solver.getPatchSpecCount("patch2", "F"), 0)
        self.assertEqual(solver.getPatchSpecCount("patch1", "A"), solver.getPatchSpecCount("patch1", "B"))
        self.assertEqual(solver.getPatchSpecCount("patch1", "A") + solver.getPatchSpecCount("patch1", "C"), 100000)
        self.assertEqual(solver.getPatchSpecCount("patch2", "D"), solver.getPatchSpecCount("patch2", "E"))
        self.assertEqual(solver.getPatchSpecCount("patch2", "D") + solver.getPatchSpecCount("patch2", "F"), 100000)

    def testWmdirect(self):
        solver = solv.Wmdirect(self.model, self.geom, self.rng)
        self._runTest(solver)
    
    def testWmrssa(self):
        solver = solv.Wmrssa(self.model, self.geom, self.rng)
        self._runTest(solver)

    def testTetexact(self):
        solver = solv.Tetexact(self.model, self.mesh, self.rng)
        self._runTest(solver)
def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(MultiVolsysTestCase))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(MultiSurfsysTestCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
