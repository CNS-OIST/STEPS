####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

import steps.geom as sgeom
from steps.geom import INDEX_DTYPE
from steps.utilities import meshio

import numpy as np

class GenPointsInTetsBugfixCase(unittest.TestCase):
    """ 
    Test whether Tetmesh::genPointsInTet and Tetmesh::genPointsInTri correctly output randomly
    selected points inside tetrahedrons. Fail if all points are the same.
    """
    def setUp(self):
        if __name__ == "__main__":
            self.mesh = meshio.importAbaqus('meshes/test.inp', 1e-7)[0]
        else:
            self.mesh = meshio.importAbaqus('genPointsInTets_bugfix_test/meshes/test.inp', 1e-7)[0]

        self.NbGen = 100

    def testGenInTets(self):
        tets = np.array([0], dtype=INDEX_DTYPE)
        counts = np.array([self.NbGen], dtype=np.uint32)
        data = np.zeros(self.NbGen * 3)

        self.mesh.genTetVisualPointsNP(tets, counts, data)

        data.shape = -1, 3

        self.assertTrue(len(set(tuple(point) for point in data)) > 1)

    def testGenInTris(self):
        tris = np.array([0], dtype=INDEX_DTYPE)
        counts = np.array([self.NbGen], dtype=np.uint32)
        data = np.zeros(self.NbGen * 3)

        self.mesh.genTriVisualPointsNP(tris, counts, data)

        data.shape = -1, 3

        self.assertTrue(len(set(tuple(point) for point in data)) > 1)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(GenPointsInTetsBugfixCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

