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

import os.path as osp

import numpy
import unittest

import mpi4py.MPI as MPI

import steps.geom as sgeom

TEST_DIR = osp.join(osp.dirname(osp.realpath(__file__)), "..", "..", "..")
MESH_DIR = osp.join(TEST_DIR, "mesh")


class LibraryTestCase(unittest.TestCase):
    def test_default_library(self):
        lib = sgeom.Library()

    def test_provide_communicator(self):
        lib = sgeom.Library(MPI.COMM_WORLD)


class TwoTetsTestCase(unittest.TestCase):
    def setUp(self):
        self.comm = MPI.COMM_WORLD
        self.comm_size = self.comm.Get_size()
        self.comm_rank = self.comm.Get_rank()
        self.library = sgeom.Library(self.comm)
        self.mesh = sgeom.DistMesh(self.library, osp.join(MESH_DIR, "2_tets.msh"))

    def test_elems_and_bounds(self):
        """check owned and total number of elems and bounds"""
        assert self.mesh is not None
        assert self.mesh.num_elems == 2 / self.comm_size
        assert self.mesh.total_num_elems == 2
        assert self.mesh.total_num_bounds == 7
        assert self.mesh.num_bounds in [3, 4]
        bounds = numpy.zeros(1)
        bounds[0] = self.mesh.num_bounds
        total_bounds = self.comm.allreduce(bounds, op=MPI.SUM)
        assert total_bounds[0] == 7


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(LibraryTestCase))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TwoTetsTestCase))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


