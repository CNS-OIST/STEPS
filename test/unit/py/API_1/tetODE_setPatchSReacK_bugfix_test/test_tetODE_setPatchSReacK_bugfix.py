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

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.solver as ssolver
import steps.rng as srng
from steps.utilities import meshio

class TetODESetPatchSReacKBugfixTestCase(unittest.TestCase):
    """ 
    Test whether _setPatchSReacK from tetODE works without throwing a segfault.
    """
    def setUp(self):
        mdl = smodel.Model()

        S1 = smodel.Spec('S1', mdl)

        ssys = smodel.Surfsys('ssys', mdl)

        smodel.SReac('SR01', ssys, slhs=[S1], srhs=[S1], kcst=1)

        if __name__ == "__main__":
            mesh = meshio.importAbaqus('meshes/test.inp', 1e-7)[0]
        else:
            mesh = meshio.importAbaqus('tetODE_setPatchSReacK_bugfix_test/meshes/test.inp', 1e-7)[0]

        comp1 = sgeom.TmComp('comp1', mesh, range(mesh.countTets()))

        patch1 = sgeom.TmPatch('patch1', mesh, mesh.getSurfTris(), comp1)
        patch1.addSurfsys('ssys')

        rng = srng.create('mt19937',512)
        rng.initialize(1234)

        self.sim = ssolver.TetODE(mdl, mesh, rng)
        self.sim.setTolerances(1.0e-3, 1.0e-3)

    def testSetPatchSReacK(self):
        self.sim.setPatchSReacK('patch1', 'SR01', 2)
        # getPatchSReacK or getTriSReacK are not implemented so we cannot check whether it was 
        # successful.


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TetODESetPatchSReacKBugfixTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
