# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2014 Okinawa Institute of Science and Technology, Japan.
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

""" Unit tests for steps.model.Spec class."""

import unittest2
import random
from steps.model import *

class SpecCreationTestCase(unittest2.TestCase):
    """ Tests for species construction. """
    def setUp(self):
        self.model = Model()
    
    def tearDown(self):
        self.model = None
    
    def testCreateNormalSpec(self):
        spec = Spec("A", self.model)
        self.assertIsInstance(spec, Spec)

    def testCreateValenceSpec(self):
        spec = Spec("A", self.model, random.randint(0,10))
        self.assertIsInstance(spec, Spec)

class SpecPropertyTestCase(unittest2.TestCase):
    """ Property tests for species. """
    def setUp(self):
        self.model = Model()
        self.spec_normal = Spec("normal", self.model)
        self.spec_valence = Spec("valence", self.model, 5)
    
    def tearDown(self):
        self.model = None
        self.spec_normal = None
        self.spec_valence = None
    
    def testSetGetSpecID(self):
        self.assertEqual(self.spec_normal.getID(), "normal")
        self.assertEqual(self.spec_valence.getID(), "valence")
        self.spec_normal.setID("new_id")
        self.assertEqual(self.spec_normal.getID(), "new_id")

    def testGetModel(self):
        self.assertEqual(self.spec_normal.getModel().this, self.model.this)

    def testSetGetValence(self):
        self.assertEqual(self.spec_valence.getValence(), 5)
        new_v = random.randint(0,10)
        self.spec_valence.setValence(new_v)
        self.assertEqual(self.spec_valence.getValence(), new_v)

def suite():
    all_tests = []
    all_tests.append(unittest2.makeSuite(SpecCreationTestCase, "test"))
    all_tests.append(unittest2.makeSuite(SpecPropertyTestCase, "test"))
    return unittest2.TestSuite(all_tests)

if __name__ == "__main__":
    unittest2.TextTestRunner(verbosity=2).run(suite())
