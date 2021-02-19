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

""" Unit tests for Membrane class and related methods."""

import unittest
import os

from steps import interface

from steps.model import *
from steps.geom import *

FILEDIR = os.path.dirname(os.path.abspath(__file__))

class membraneTests(unittest.TestCase):
    """Test Membrane class."""
    def setUp(self):
        mdl = Model()
        with mdl:
            vsys1, vsys2 = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

        self.mesh = TetMesh.Load(os.path.join(FILEDIR, 'meshes', 'cyl_len10_diam1'))

        with self.mesh:
            center1 = self.mesh.tets[0, 0, self.mesh.bbox.max.z / 2]
            self.c1tets = TetList([tet for tet in self.mesh.tets if tet.center.z > 0 and tet.idx != center1.idx])
            c2tets = self.mesh.tets - self.c1tets - TetList([center1])
            comp1 = Compartment.Create(self.c1tets, mdl.vsys1)
            comp1bis = Compartment.Create([center1], mdl.vsys1)
            comp2 = Compartment.Create(c2tets, mdl.vsys2)
            ptris = self.c1tets.surface & c2tets.surface
            self.ptris = TriList(tri for tri in ptris if tri.center.x > 0)
            self.ptris2 = ptris - self.ptris
            patch = Patch.Create(self.ptris, comp1, comp2, mdl.ssys)
            patch2 = Patch.Create(self.ptris2, comp1, comp2, mdl.ssys)

            membrane = Membrane.Create([patch, patch2], opt_method = 1)

        self.mdl = mdl
        self.memb = membrane


    def testPreventSystemAdding(self):
        """Check that surface systems cannot be added to membranes."""

        with self.assertRaises(NotImplementedError):
            self.memb.addSystem(self.mdl.ssys)

    def testProperties(self):
        """Test accessing membrane properties."""

        smemb = self.memb._getStepsObjects()[0]
        self.assertEqual(self.memb.name, smemb.getID())
        self.assertEqual(self.memb.open, smemb.open())
        self.assertAlmostEqual(self.memb.Area, self.ptris.Area + self.ptris2.Area)





def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(membraneTests, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


