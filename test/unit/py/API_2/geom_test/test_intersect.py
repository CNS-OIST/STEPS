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

""" Unit test to validate intersect method """

import unittest
import os.path as osp
import numpy as np

from steps import interface
from steps.geom import *
from steps.rng import *

TEST_DIR = osp.join(osp.dirname(osp.realpath(__file__)), "..", "..", "..", "..")
MESH_DIR = osp.join(TEST_DIR, "mesh")


class IntersectTests(unittest.TestCase):
    def setUp(self):
        """ Setup meshes -> STEPS3/4 """
        # cube.msh is a cube of edge length 1 and edge corner at [0,0,0]
        self.distMesh = DistMesh(osp.join(MESH_DIR, "cube.msh"))
        self.tetMesh = TetMesh.LoadGmsh(osp.join(MESH_DIR, "cube.msh"))

    def testIntersect_compareSTEPS4WithSTEPS3(self):
        """ comparison between STEPS3 & STEPS4 and other sanity checks """

        # Test 1 : a line with 1 segment
        # Create a line from [0,0,0] to [1,1,1]
        # From one corner of the cube to the other
        pts = np.array([[0,0,0],[1,1,1]], dtype=float, order='C')

        # For every segment we pass, we get a vector/list of pairs of (tet, intersection fraction)
        intersect_distMesh = self.distMesh.stepsMesh.intersect(pts)
        intersect_tetMesh = self.tetMesh.stepsMesh.intersect(pts)
        self.check_intersect(intersect_distMesh, intersect_tetMesh, 1)

        # Test 2 : a line with 2 segments
        pts = np.array([[0,0,0],[0.5,0.5,0.5],[1,1,1]], dtype=float, order='C')
        intersect_distMesh = self.distMesh.stepsMesh.intersect(pts)
        intersect_tetMesh = self.tetMesh.stepsMesh.intersect(pts)
        self.check_intersect(intersect_distMesh, intersect_tetMesh, 2)

    def check_intersect(self, intersect_distMesh, intersect_tetMesh, segments):
        self.assertEqual(len(intersect_distMesh), segments)
        self.assertEqual(len(intersect_distMesh), len(intersect_tetMesh))
        self.assertListEqual(intersect_distMesh, intersect_tetMesh)
        not_empty = 0
        for segment in intersect_distMesh:
            for p in segment:
                # p[0] : tet ID
                # p[1] : intersection fraction
                self.assertLess(p[0], self.distMesh.stepsMesh.total_num_elems)
                self.assertLessEqual(p[1], 1)
                self.assertGreaterEqual(p[0], 0)
                self.assertGreaterEqual(p[1], 0)
                not_empty += 1
        self.assertGreater(not_empty, 0, 
            "No pairs of (tet, fract) returned, even if the line is inside the mesh.")


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(IntersectTests, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
