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
        with self.distMesh.asLocal():
            eps = 1e-3

            # Test 1 : a line with 1 segment -> from one corner of the cube to the other
            # The eps is needed because the points of the segments are checked in which tetrahedron
            # they belong to, and for this reason we want to avoid having them on the corners/faces.
            pts = np.array([[0+eps,0+eps,0+eps],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            # STEPS4: Start point outside the mesh
            pts_start_out = np.array([[-1,-1,-1],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            # STEPS4: End point outside the mesh
            pts_end_out = np.array([[0+eps,0+eps,0+eps],[2,2,2]], dtype=float, order='C')
            # STEPS4: Start/End points outside the mesh
            pts_both_out = np.array([[-1,-1,-1],[2,2,2]], dtype=float, order='C')

            # For every segment we pass, we get a vector/list of pairs of (tet, intersection fraction)
            intersect_distMesh = self.distMesh.intersect(pts)
            intersect_tetMesh = self.tetMesh.intersect(pts)
            self.check_intersect(intersect_distMesh, intersect_tetMesh, 1)

            # Test STEPS4, when start/end points are outside the mesh
            # check that they are passing through the same tets
            intersect_distMesh_s = self.distMesh.intersect(pts_start_out)
            self.check_crossing_tets(intersect_distMesh_s, intersect_distMesh)
            intersect_distMesh_e = self.distMesh.intersect(pts_end_out)
            self.check_crossing_tets(intersect_distMesh_e, intersect_distMesh)
            intersect_distMesh_b = self.distMesh.intersect(pts_both_out)
            self.check_crossing_tets(intersect_distMesh_b, intersect_distMesh)

            # Test 2 : a line with 2 segments
            pts = np.array([[0+eps,0+eps,0+eps],[0.5,0.5,0.5],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            pts_start_out = np.array([[-1,-1,-1],[0.5,0.5,0.5],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            pts_end_out = np.array([[0+eps,0+eps,0+eps],[0.5,0.5,0.5],[2,2,2]], dtype=float, order='C')
            pts_both_out = np.array([[-1,-1,-1],[0.5,0.5,0.5],[2,2,2]], dtype=float, order='C')

            intersect_distMesh = self.distMesh.intersect(pts)
            intersect_tetMesh = self.tetMesh.intersect(pts)
            self.check_intersect(intersect_distMesh, intersect_tetMesh, 2)

            intersect_distMesh_s = self.distMesh.intersect(pts_start_out)
            self.check_crossing_tets(intersect_distMesh_s, intersect_distMesh)
            intersect_distMesh_e = self.distMesh.intersect(pts_end_out)
            self.check_crossing_tets(intersect_distMesh_e, intersect_distMesh)
            intersect_distMesh_b = self.distMesh.intersect(pts_both_out)
            self.check_crossing_tets(intersect_distMesh_b, intersect_distMesh)

    def testIntersect_independentSegments(self):
        """ TODO """
        with self.distMesh.asLocal():
            eps = 1e-3

            # Test 1 : a line with 1 segment -> from one corner of the cube to the other
            # The eps is needed because the points of the segments are checked in which tetrahedron
            # they belong to, and for this reason we want to avoid having them on the corners/faces.
            pts = np.array([[0.1e-6, 0.1e-6, 0.1e-6],[0.9e-6, 0.9e-6, 0.9e-6]], dtype=float, order='C')

            # For every segment we pass, we get a vector/list of pairs of (tet, intersection fraction)
            intersect_distMesh = self.distMesh.intersectIndependentSegments(pts)

    def check_intersect(self, intersect_distMesh, intersect_tetMesh, segments):
        self.assertEqual(len(intersect_distMesh), segments)
        self.assertEqual(len(intersect_distMesh), len(intersect_tetMesh))

        dist_set = [sorted([(tet.idx, round(rat, 3)) for tet, rat in seg]) for seg in intersect_distMesh]
        tet_set = [sorted([(tet.idx, round(rat, 3)) for tet, rat in seg]) for seg in intersect_tetMesh]

        self.assertEqual(dist_set, tet_set)
        not_empty = 0
        for segment in intersect_distMesh:
            for p in segment:
                # p[0] : tet ID
                # p[1] : intersection fraction
                self.assertLess(p[0].idx, len(self.distMesh.tets))
                self.assertLessEqual(p[1], 1)
                self.assertGreaterEqual(p[0].idx, 0)
                self.assertGreaterEqual(p[1], 0)
                not_empty += 1
        self.assertGreater(not_empty, 0, 
            "No pairs of (tet, fract) returned, even if the line is inside the mesh.")
    
    def check_crossing_tets(self, intersect1, intersect2):
        tets1 = sorted([tet.idx for seg in intersect1 for tet,_ in seg])
        tets2 = sorted([tet.idx for seg in intersect2 for tet,_ in seg])
        self.assertEqual(tets1, tets2)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(IntersectTests, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
