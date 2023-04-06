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

# Test utility functions in steps.utilities.meshctrl.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.geom as sgeom
from steps.utilities import meshio
from steps.utilities import meshctrl

import numpy as np

class MeshctrlUtilityTestCase(unittest.TestCase):
    """ 
    Test utility functions in steps.utilities.meshctrl.
    """
    def setUp(self):
        if __name__ == "__main__":
            self.left_mesh = meshio.importAbaqus('meshes/comp1.inp', 1e-6)[0]
            self.right_mesh = meshio.importAbaqus('meshes/comp2.inp', 1e-6)[0]
            self.combine_mesh = meshio.importAbaqus('meshes/2comps.inp', 1e-6)[0]
        else:
            self.left_mesh = meshio.importAbaqus('utilities_meshctrl_test/meshes/comp1.inp', 1e-6)[0]
            self.right_mesh = meshio.importAbaqus('utilities_meshctrl_test/meshes/comp2.inp', 1e-6)[0]
            self.combine_mesh = meshio.importAbaqus('utilities_meshctrl_test/meshes/2comps.inp', 1e-6)[0]
        self.test_comp = sgeom.TmComp("test_comp", self.combine_mesh, range(12))

    def testFindOverlapTris(self):
        overlap_tris = meshctrl.findOverlapTris(self.combine_mesh, range(12), range(12, 24))

        self.assertEqual(len(overlap_tris), 2)
        for tri in overlap_tris:
            vertices = self.combine_mesh.getTri(tri)
            for v in vertices:
                coords = self.combine_mesh.getVertex(v)
                self.assertTrue(np.isclose(coords[0], 0.0, rtol=1e-12, atol=1e-15))

    def testFindOverlapSurfTris(self):
        overlap_couplings = meshctrl.findOverlapSurfTris(self.left_mesh, self.right_mesh)

        self.assertEqual(len(overlap_couplings), 2)
        for coupling_data in overlap_couplings:
            surftri_1 = coupling_data[1]
            surftri_2 = coupling_data[3]
            vertices_1 = self.left_mesh.getTri(surftri_1)
            vertices_2 = self.left_mesh.getTri(surftri_2)
            coords_1 = np.sort(np.array([self.left_mesh.getVertex(v) for v in vertices_1]))
            coords_2 = np.sort(np.array([self.right_mesh.getVertex(v) for v in vertices_2]))
            print(coords_1)
            print(coords_2)
            self.assertTrue(np.allclose(coords_1, coords_2))

    def testFindSurfTris(self):
        comp_surf_tris = meshctrl.findSurfTrisInComp(self.combine_mesh, self.test_comp)
        tet_surf_tris = meshctrl.findSurfTrisInTets(self.combine_mesh, range(12))
        self.assertTrue(np.array_equal(np.sort(comp_surf_tris), np.sort(comp_surf_tris)))
        self.assertEqual(len(tet_surf_tris), 10)

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(MeshctrlUtilityTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

