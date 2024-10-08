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

""" Unit test to validate intersect method [distributed mesh] """

import unittest
import os.path as osp
import numpy as np
import mpi4py

from steps import interface
from steps.geom import *
from steps.rng import *

TEST_DIR = osp.join(osp.dirname(osp.realpath(__file__)), "..", "..", "..", "..")
MESH_DIR = osp.join(TEST_DIR, "mesh")


class distIntersectTests(unittest.TestCase):
    def setUp(self):
        """ Setup meshes -> STEPS3/4 """
        # cube.msh is a cube of edge length 1 and edge corner at [0,0,0]
        self.tetMesh = TetMesh.LoadGmsh(osp.join(MESH_DIR, "cube.msh"))

    def testIntersectDistributedLocal_n2(self):
        """ Test STEPS 4 intersect on distributed mesh """
        # splitMesh = DistMesh(osp.join(MESH_DIR, "cube_split_2/cube"))
        splitMesh = DistMesh(osp.join(MESH_DIR, "3_tets.msh"))
        MPI_rank = mpi4py.MPI.COMM_WORLD.Get_rank()
        MPI_size = mpi4py.MPI.COMM_WORLD.Get_size()
        
        with splitMesh.asLocal():
            pts = np.array([[0.2,0.1, 0.1],[0.1, 0.1, 0.1]], dtype=float, order='C')
            ans = self.sort_and_round_intersections(splitMesh.intersect(pts))
            self.assertEqual(ans, [[(0, 1.0)]] if MPI_rank == 1 else [[]])

            pts = np.array([[0.2,0.1, 0.1],[-0.2, 0.1, 0.1]], dtype=float, order='C')
            ans = self.sort_and_round_intersections(splitMesh.intersect(pts))
            self.assertEqual(ans, [[(0, 0.5), (2, 0.5)]] if MPI_rank == 1 else [[]])

            pts = np.array([[0.2,0.1, -0.1],[0.2, 0.1, 0.3]], dtype=float, order='C')
            ans = self.sort_and_round_intersections(splitMesh.intersect(pts))
            self.assertEqual(ans, [[(0, 0.75)]] if MPI_rank == 1 else [[(1, 0.25)]])

            pts = np.array([[0.2,0.1, -0.1],[0.1, 0.1, 0]], dtype=float, order='C')
            ans = self.sort_and_round_intersections(splitMesh.intersect(pts))
            self.assertEqual(ans, [[]] if MPI_rank == 1 else [[(1, 1.0)]])
            
            pts = np.array([[0.2,0.1, -0.1],[0, 0.1, 0.1]], dtype=float, order='C')
            ans = self.sort_and_round_intersections(splitMesh.intersect(pts))
            self.assertEqual(ans, [[(0, 0.5)]] if MPI_rank == 1 else [[(1, 0.5)]])

            pts = np.array([[0.2,0.1, -0.1],[-0.2, 0.1, 0.3]], dtype=float, order='C')
            ans = self.sort_and_round_intersections(splitMesh.intersect(pts))
            self.assertEqual(ans, [[(0, 0.25), (2, 0.5)]] if MPI_rank == 1 else [[(1, 0.25)]])

    def sort_and_round_intersections(self, intersec):
        """ Put the intersections in a standard form so they can be easily compared """
        return [sorted([(tet.toGlobal().idx, round(rat, 3)) for tet, rat in seg]) for seg in intersec]

    def testIntersectDistributedSTEPS3vsSTEPS4_n2(self):
        """ Comparison between STEPS3 & distributed STEPS4 and other sanity checks """

        splitMesh = DistMesh(osp.join(MESH_DIR, "cube_split_2/cube"))
        with splitMesh.asLocal():
            eps = 1e-3

            # Test 1 : a line with 1 segment -> from one corner of the cube to the other
            # The eps is needed because the points of the segments are checked in which tetrahedron
            # they belong to, and for this reason we want to avoid having them on the corners/faces.
            pts = np.array([[0+eps,0+eps,0+eps],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            # STEPS4: Start point outside the mesh
            pts_start_out = np.array([[-1,-1,-1],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            # STEPS4: End point outside the mesh
            pts_end_out = np.array([[0+eps,0+eps,0+eps],[2-eps,2-eps,2-eps]], dtype=float, order='C')
            # STEPS4: Start/End points outside the mesh
            pts_both_out = np.array([[-1+eps,-1+eps,-1+eps],[2,2,2]], dtype=float, order='C')

            # For every segment we pass, we get a vector/list of pairs of (tet, intersection fraction)
            intersect_tetMesh = self.tetMesh.intersect(pts)
            intersect_distMesh = splitMesh.intersect(pts)
            self.check_tets_and_ratios(splitMesh, intersect_distMesh, intersect_tetMesh)

            # Test raw and local kwyword arguments
            self.assertTrue(all(all(tet.isLocal() for tet, rat in seg) for seg in intersect_distMesh))

            raw_intersect_distMesh = splitMesh.intersect(pts, raw=True)
            self.assertEqual(raw_intersect_distMesh, [[(tet.idx, rat) for tet, rat in seg] for seg in intersect_distMesh])

            global_intersect_distMesh = splitMesh.intersect(pts, local=False)
            self.assertTrue(all(all(not tet.isLocal() for tet, rat in seg) for seg in global_intersect_distMesh))
            self.assertEqual(global_intersect_distMesh, [[(tet.toGlobal(), rat) for tet, rat in seg] for seg in intersect_distMesh])

            # Test STEPS4, when start/end points are outside the mesh
            # check that they are passing through the same tets
            intersect_distMesh = splitMesh.intersect(pts_start_out)
            self.check_crossing_tets(splitMesh, intersect_distMesh, intersect_tetMesh)
            intersect_distMesh = splitMesh.intersect(pts_end_out)
            self.check_crossing_tets(splitMesh, intersect_distMesh, intersect_tetMesh)
            intersect_distMesh = splitMesh.intersect(pts_both_out)
            self.check_crossing_tets(splitMesh, intersect_distMesh, intersect_tetMesh)

            # Test 2 : a line with 2 segments
            pts = np.array([[0+eps,0+eps,0+eps],[0.5,0.5,0.5],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            pts_start_out = np.array([[-1+eps,-1,-1+eps],[0.5,0.5,0.5],[1-eps,1-eps,1-eps]], dtype=float, order='C')
            pts_end_out = np.array([[0+eps,0+eps,0+eps],[0.5,0.5,0.5],[2-eps,2-eps,2-eps]], dtype=float, order='C')
            pts_both_out = np.array([[-1+eps,-1,-1+eps],[0.5,0.5,0.5],[2,2,2]], dtype=float, order='C')

            intersect_tetMesh = self.tetMesh.intersect(pts)
            intersect_distMesh = splitMesh.intersect(pts)
            self.check_tets_and_ratios(splitMesh, intersect_distMesh, intersect_tetMesh)

            intersect_distMesh = splitMesh.intersect(pts_start_out)
            self.check_crossing_tets(splitMesh, intersect_distMesh, intersect_tetMesh)
            intersect_distMesh = splitMesh.intersect(pts_end_out)
            self.check_crossing_tets(splitMesh, intersect_distMesh, intersect_tetMesh)
            intersect_distMesh = splitMesh.intersect(pts_both_out)
            self.check_crossing_tets(splitMesh, intersect_distMesh, intersect_tetMesh)

    def check_tets_and_ratios(self, splitMesh, intersect_distMesh, intersect_tetMesh):
        # number of segments
        nsegs = 0
        for seg in intersect_tetMesh:
            nsegs += 1
        
        # The choice of barycenters is because the split mesh has different tet ids compared to
        # the non-split mesh. Therefore, the comparison is done through a common property.
        
        # non-split mesh: gather tet barycenters and ratios
        tet_list = [(tuple(tet.center), round(ratio,4)) for seg in intersect_tetMesh for tet,ratio in seg]

        # split mesh: gather tet barycenters and ratios
        map = {tet.idx: tuple(tet.center) for tet in splitMesh.tets}
        dist_list = [(map[tet.idx], round(ratio,4)) for seg in intersect_distMesh for tet,ratio in seg]
        # gather tuples across ranks and flatten the list [after allgather, same container across ranks]
        dist_list = mpi4py.MPI.COMM_WORLD.allgather(dist_list)
        dist_list = [tuple for task in dist_list for tuple in task]

        self.assertEqual(sorted(dist_list), sorted(tet_list))

        # Check that ratios per segment sum to 1
        ratio_tetMesh = [ratio for seg in intersect_tetMesh for _,ratio in seg]
        assert np.allclose(sum(ratio_tetMesh), 1.0*nsegs, rtol=1e-4, atol=1e-4)

        ratio_distMesh = [pair[1] for pair in dist_list]
        assert np.allclose(sum(ratio_distMesh), 1.0*nsegs, rtol=1e-4, atol=1e-4)

    def check_crossing_tets(self, splitMesh, intersect_distMesh, intersect_tetMesh):
        # non-split mesh: gather tet barycenters
        tet_sorted_set = sorted(list(set([tuple(tet.center) for seg in intersect_tetMesh for tet,_ in seg])))

        # split mesh: gather tet barycenters
        map = {tet.idx: tuple(tet.center) for tet in splitMesh.tets}
        dist_sorted_set= sorted(list(set([map[tet.idx] for seg in intersect_distMesh for tet,_ in seg])))
        # gather across ranks and flatten [after allgather, same container across ranks]
        dist_sorted_set = mpi4py.MPI.COMM_WORLD.allgather(dist_sorted_set)
        dist_sorted_set = sorted(list(set(([tuple for task in dist_sorted_set for tuple in task]))))
        
        self.assertEqual(tet_sorted_set, dist_sorted_set)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(distIntersectTests, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
