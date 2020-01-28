# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import print_function
import unittest2
from steps.utilities import meshio

class GmshImportTest(unittest2.TestCase):
    """ Test if a Gmsh mesh with multiple tags can be imported correctly. """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testImportGmsh2ASCII(self):
        if __name__ == "__main__":
            mesh_file = "meshes/two_tag_sphere_v2ascII.msh"
        else:
            mesh_file = "gmsh_multitag_test/meshes/two_tag_sphere_v2ascII.msh"
        mesh, node_proxy, tet_proxy, tri_proxy = meshio.importGmsh(mesh_file, 1)
        tet_groups = tet_proxy.getGroups()
        self.assertIn((0, 'inner'), tet_groups.keys())
        self.assertIn((0, 'outer'), tet_groups.keys())
        self.assertIn((0, 1), tet_groups.keys())
        self.assertIn((0, 2), tet_groups.keys())
        self.assertIn((1, 3), tet_groups.keys())
        self.assertIn((1, 4), tet_groups.keys())

    def testImportGmsh4ASCIIFailCase(self):
        if __name__ == "__main__":
            mesh_file = "meshes/two_tag_sphere_v4ascII.msh"
        else:
            mesh_file = "gmsh_multitag_test/meshes/two_tag_sphere_v4ascII.msh"
        self.assertRaises(Exception, meshio.importGmsh, mesh_file, 1)

    def testImportGmshBinaryFailCase(self):
        if __name__ == "__main__":
            mesh_file = "meshes/two_tag_sphere_v4binary.msh"
        else:
            mesh_file = "gmsh_multitag_test/meshes/two_tag_sphere_v4binary.msh"
        self.assertRaises(Exception, meshio.importGmsh, mesh_file, 1)
        

def suite():
    all_tests = []
    all_tests.append(unittest2.makeSuite(GmshImportTest, "test"))
    return unittest2.TestSuite(all_tests)

if __name__ == "__main__":
    unittest2.TextTestRunner(verbosity=2).run(suite())