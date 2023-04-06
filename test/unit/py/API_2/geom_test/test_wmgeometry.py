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

""" Unit tests for Geometry class and related methods."""

import unittest
import os

from steps import interface

from steps.geom import *

class wmGeometryTests(unittest.TestCase):
    """Test well-mixed geometry classes."""
    def setUp(self):
        self.geom = Geometry()

        with self.geom:
            comp1 = Compartment.Create(None, 1)
            comp2 = Compartment.Create(None, 2)
            comp3 = Compartment.Create(None, 3)

            patch1 = Patch.Create(comp1, comp2, None, 1)
            patch2 = Patch.Create(comp2, comp3, None, 2)
            patch3 = Patch.Create(comp3, comp1, None, 3)

    def testFromSteps(self):
        """Test loading of Geometry from a steps object."""

        geom2 = Geometry._FromStepsObject(self.geom._getStepsObjects()[0])

        comps = {c.name: c for c in geom2.ALL(Compartment)}
        patches = {c.name: c for c in geom2.ALL(Patch)}

        self.assertEqual(len(comps), 3)
        self.assertEqual(len(patches), 3)

        for c in self.geom.ALL(Compartment):
            c2 = comps[c.name]
            self.assertEqual(c2.Vol, c.Vol)
            self.assertEqual(set(c2.sysNames), set(c.sysNames))
            self.assertEqual(c2._getStepsObjects()[0].getID(), c._getStepsObjects()[0].getID())

        for p in self.geom.ALL(Patch):
            p2 = patches[p.name]
            self.assertEqual(p2.Area, p.Area)
            self.assertEqual(set(p2.sysNames), set(p.sysNames))
            self.assertEqual(p2.innerComp.name, p.innerComp.name)
            self.assertEqual(p2.outerComp.name, p.outerComp.name)





def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(wmGeometryTests, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

