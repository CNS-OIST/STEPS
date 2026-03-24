####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

# Test utility functions in steps.utilities.geom_decompose.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest
import os

from steps.API_2.geom import TetMesh, Patch, Compartment, _getTriPartitionFromTet
from steps.API_2.model import Model, Species, VolumeSystem, Diffusion


FILEDIR = os.path.dirname(os.path.abspath(__file__))

class GeomDecomposeTestCase(unittest.TestCase):
    """
    Test utility functions in steps.utilities.geom_decompose.
    """
    def setUp(self):
        self.mesh = TetMesh.LoadGmsh(os.path.join(FILEDIR, '../../../../mesh/box.msh'), scale=1e-6)
        self.model = Model()
        with self.model:
            vsys = VolumeSystem.Create()
            SA = Species.Create()
            with vsys:
                Diffusion(SA, 0)
            self.vsys = vsys


    def testPartitionTris(self):
        mesh = self.mesh
        with mesh:
            # Create one compartment per tet
            for tet in mesh.tets:
                Compartment(tet.toList(), self.vsys)
            # Create one patch per triangle between two compartments
            for tri in mesh.tris - mesh.surface:
                comps = [tet.comp for tet in tri.tetNeighbs]
                Patch(tri.toList(), *comps)

         # One host per tetrahedron
        tet_hosts = [i for i in range(len(mesh.tets))]

        # Partition triangles
        tri_hosts = _getTriPartitionFromTet(mesh, tet_hosts)

        # Assert that tet_hosts are all 0
        self.assertTrue(all(host == 0 for host in tet_hosts))

        # Assert that tri_hosts are all 0
        self.assertTrue(all(host == 0 for host in tri_hosts.values()))


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(GeomDecomposeTestCase))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

