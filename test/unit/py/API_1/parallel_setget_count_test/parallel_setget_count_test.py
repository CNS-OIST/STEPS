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

class ParallelSetGetCountCase(unittest.TestCase):
    """ 
    Test if set/getTet/TriCount in TetOpSplit works as expected.
    Used to detect bug in https://github.com/CNS-OIST/HBP_STEPS/issues/229
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
            self.mesh = meshio.importAbaqus('parallel_diff_sel_test/meshes/test_mesh.inp', 1e-7)[0]

        self.tmcomp = sgeom.TmComp('comp', self.mesh, range(self.mesh.ntets))
        self.tmcomp.addVolsys('vsys')
        self.surf_tris = self.mesh.getSurfTris()
        self.tmpatch = sgeom.TmPatch('patch', self.mesh, self.surf_tris, icomp = self.tmcomp)
        self.tmpatch.addSurfsys('ssys')

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)
        
    def tearDown(self):
        self.model = None
        self.mesh = None
        self.rng = None

    def testSetGetTetCount(self):
        tet_hosts = gd.binTetsByAxis(self.mesh, steps.mpi.nhosts)
        tri_hosts = gd.partitionTris(self.mesh, tet_hosts, self.surf_tris)
        solver = solv.TetOpSplit(self.model, self.mesh, self.rng, solv.EF_NONE, tet_hosts, tri_hosts)
        for tet in range(10):
            solver.setTetSpecCount(tet, 'A', tet)
            get_count = solver.getTetSpecCount(tet, 'A')
            self.assertEqual(get_count, tet)

    def testSetGetTriCount(self):
        tet_hosts = gd.binTetsByAxis(self.mesh, steps.mpi.nhosts)
        tri_hosts = gd.partitionTris(self.mesh, tet_hosts, self.surf_tris)
        solver = solv.TetOpSplit(self.model, self.mesh, self.rng, solv.EF_NONE, tet_hosts, tri_hosts)
        for tri in self.surf_tris:
            solver.setTriSpecCount(tri, 'A', tri)
            get_count = solver.getTriSpecCount(tri, 'A')
            self.assertEqual(get_count, tri)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(ParallelSetGetCountCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
