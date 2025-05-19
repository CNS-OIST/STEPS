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

""" Unit tests for tetmesh class and related methods."""

import unittest

from steps import interface

from steps.geom import *

import os

FILEDIR = os.path.dirname(os.path.abspath(__file__))

class BboxDistMeshTest(unittest.TestCase):
    """Test whether the bounding box is computed correctly with split meshes."""
    def testBoundingBox(self):
        if DistMesh._use_gmsh():
            mesh = DistMesh(os.path.join(FILEDIR, 'meshes', 'split_cube'))

            positions = [
                [0.7210276 , 0.76670546, 0.96564512],
                [0.89970559, 0.27098829, 0.93006672],
                [0.51967092, 0.90021152, 0.5356602 ],
                [0.03331348, 0.45842372, 0.18638614],
                [0.49220902, 0.58731909, 0.51503226],
            ]

            # Should not raise any exceptions
            for pos in positions:
                mesh.tets[pos]

            with self.assertRaises(KeyError):
                mesh.tets[2, 0, 0]

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(BboxDistMeshTest))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
