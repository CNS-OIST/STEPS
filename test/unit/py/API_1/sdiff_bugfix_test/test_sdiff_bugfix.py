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

# Unit test for bugfix https://github.com/CNS-OIST/HBP_STEPS/pull/144

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

from steps.geom import UNKNOWN_TRI
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solv
from steps.utilities import meshio

class SDiffBndSpeciesIDTestCase(unittest.TestCase):
    """ Test if the species ids of a molecule on two patches with different vsys between a surface boundary is mapped correctly during a diffusion event. """
    def setUp(self):
        self.DCST = 0.08e-12

        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)
        X = smodel.Spec('X', self.model)

        self.ssys1 = smodel.Surfsys('ssys1', self.model)
        self.ssys2 = smodel.Surfsys('ssys2', self.model)
    
        self.sreac = smodel.SReac('sreac', self.ssys1, slhs = [A], srhs = [A],  kcst = 1e5)
        self.diff1 = smodel.Diff('diffX1', self.ssys1, X, self.DCST)
        self.diff2 = smodel.Diff('diffX2', self.ssys2, X, self.DCST)

        if __name__ == "__main__":
            self.mesh = meshio.loadMesh('meshes/coin_10r_1h_13861')[0]
        else:
            self.mesh = meshio.loadMesh('sdiff_bugfix_test/meshes/coin_10r_1h_13861')[0]

        ntets = self.mesh.countTets()
        self.comp = sgeom.TmComp('cyto', self.mesh, range(ntets))

        alltris = self.mesh.getSurfTris()

        patchA_tris = []
        patchB_tris = []
        patchA_bars = set()
        patchB_bars = set()

        for t in alltris:
            vert0, vert1, vert2 = self.mesh.getTri(t)
            if (self.mesh.getVertex(vert0)[2] > 0.0 \
                and self.mesh.getVertex(vert1)[2] > 0.0 \
                and self.mesh.getVertex(vert2)[2] > 0.0):
                if self.mesh.getTriBarycenter(t)[0] > 0.0:
                    patchA_tris.append(t)
                    bar = self.mesh.getTriBars(t)
                    patchA_bars.add(bar[0])
                    patchA_bars.add(bar[1])
                    patchA_bars.add(bar[2])
                else:
                    patchB_tris.append(t)
                    bar = self.mesh.getTriBars(t)
                    patchB_bars.add(bar[0])
                    patchB_bars.add(bar[1])
                    patchB_bars.add(bar[2])

        self.patchA = sgeom.TmPatch('patchA', self.mesh, patchA_tris, icomp = self.comp)
        self.patchA.addSurfsys('ssys1')
        self.patchB = sgeom.TmPatch('patchB', self.mesh, patchB_tris, icomp = self.comp)
        self.patchB.addSurfsys('ssys2')

        # Find the set of bars that connect the two patches as the intersecting bars
        barsDB = patchA_bars.intersection(patchB_bars)
        barsDB=list(barsDB)

        # Create the surface diffusion boundary
        self.diffb = sgeom.SDiffBoundary('sdiffb', self.mesh, barsDB, [self.patchA, self.patchB])

        ctetidx = self.mesh.findTetByPoint([0.0, 0.0, 0.5e-6])
        ctet_trineighbs = self.mesh.getTetTriNeighb(ctetidx)
        self.ctri_idx = UNKNOWN_TRI
        for t in ctet_trineighbs: 
            if t in patchA_tris + patchB_tris:
                self.ctri_idx = t

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
        self.solver = solv.Tetexact(self.model, self.mesh, self.rng)
    
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None
        self.solver = None
    
    def testCrossBoundarySDiffEvent(self):
        self.solver.setTriSpecCount(self.ctri_idx, 'X', 1000)
        self.solver.setSDiffBoundarySpecDiffusionActive('sdiffb', 'X', True)
        self.solver.setSDiffBoundarySpecDcst('sdiffb', 'X', 0.008e-12 , 'patchA')
        self.solver.run(1)
        self.assertEqual(self.solver.getPatchSpecCount("patchA", "A"), 0)
        self.assertEqual(self.solver.getPatchSpecCount("patchA", "X") + self.solver.getPatchSpecCount("patchB", "X"), 1000)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(SDiffBndSpeciesIDTestCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
