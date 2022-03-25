####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

class DiffSelTestCase(unittest.TestCase):
    """ 
    Test if diffusion selector works as expected.
    Used to detect bug in https://github.com/CNS-OIST/HBP_STEPS/issues/157
    """
    def setUp(self):
        DCST = 0.08e-10
        self.model = smodel.Model()
        A = smodel.Spec('A', self.model)

        self.vsys = smodel.Volsys('vsys', self.model)

        self.diff = smodel.Diff("diff", self.vsys, A, DCST)
    
        if __name__ == "__main__":
            self.mesh = meshio.importAbaqus('meshes/test_mesh.inp', 1e-7)[0]
        else:
            self.mesh = meshio.importAbaqus('parallel_diff_sel_test/meshes/test_mesh.inp', 1e-7)[0]

        self.tmcomp = sgeom.TmComp('comp', self.mesh, range(self.mesh.ntets))
        self.tmcomp.addVolsys('vsys')

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None

    def testDiffSel(self):
        tet_hosts = gd.binTetsByAxis(self.mesh, steps.mpi.nhosts)
        solver = solv.TetOpSplit(self.model, self.mesh, self.rng, solv.EF_NONE, tet_hosts)
        solver.setCompCount('comp', 'A', 1)
        if __name__ == "__main__":
            print("Running...")
        solver.run(0.001)
        if __name__ == "__main__":
            print("")
            print("A:", solver.getCompCount("comp", "A"))
            print("steps:", solver.getNSteps())
        self.assertEqual(solver.getCompCount("comp", "A"), 1)
        self.assertNotEqual(solver.getNSteps(), 0)

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DiffSelTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
