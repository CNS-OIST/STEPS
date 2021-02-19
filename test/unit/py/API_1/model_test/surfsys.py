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

""" Unit tests for steps.model.Surfsys class."""

import unittest
import random
from steps.model import *

class SurfsysCreationTestCase(unittest.TestCase):
    """ Tests for Surfsys construction. """
    def setUp(self):
        self.model = Model()
    
    def tearDown(self):
        self.model = None
    
    def testCreateSurfsys(self):
        ssys = Surfsys("ssys", self.model)
        self.assertIsInstance(ssys, Surfsys)


class SurfsysPropertyTestCase(unittest.TestCase):
    """ Property tests for Surfsys. """
    def setUp(self):
        self.model = Model()
        self.ssys = Surfsys("ssys", self.model)
    
    def tearDown(self):
        self.model = None
        self.ssys = None
    
    def testSetGetID(self):
        self.assertEqual(self.ssys.getID(), "ssys")
        self.ssys.setID("new_id")
        self.assertEqual(self.ssys.getID(), "new_id")

    def testGetModel(self):
        self.assertEqual(self.ssys.getModel().this, self.model.this)

class SurfsysSReacTestCase(unittest.TestCase):
    """ SReac tests for Surfsys. """
    def setUp(self):
        self.model = Model()
        self.ssys = Surfsys("ssys", self.model)
        self.spec_a = Spec("A", self.model)
        self.spec_b = Spec("B", self.model)
        self.spec_c = Spec("C", self.model)
        self.sreac1 = SReac("sreac1", self.ssys, ilhs = [self.spec_a], slhs = [self.spec_b], srhs = [self.spec_c], kcst = 1e-6)
        self.sreac2 = SReac("sreac2", self.ssys, ilhs = [self.spec_a], slhs = [self.spec_b], srhs = [self.spec_c], orhs = [self.spec_a], kcst = 1e-6)
    
    def tearDown(self):
        self.model = None
        self.ssys = None
        self.spec_a = None
        self.spec_b = None
        self.spec_c = None
        self.sreac1 = None
        self.sreac2 = None
    
    def testGetSReac(self):
        self.assertEqual(self.ssys.getSReac("sreac1").this, self.sreac1.this)
    
    def testDelSReac(self):
        self.ssys.delSReac("sreac1")
        self.assertRaises(Exception, self.ssys.getSReac, "sreac1")

    def testGetAllSReacs(self):
        stored_sreacs = [sreac.this for sreac in self.ssys.getAllSReacs()]
        expected_sreacs = [self.sreac1.this, self.sreac2.this]
        self.assertCountEqual(stored_sreacs, expected_sreacs)

class SurfsysDiffTestCase(unittest.TestCase):
    """ Diff tests for Surfsys. """
    def setUp(self):
        self.model = Model()
        self.ssys = Surfsys("ssys", self.model)
        self.spec_a = Spec("A", self.model)
        self.spec_b = Spec("B", self.model)
        self.diff1 = Diff("diff1", self.ssys, self.spec_a, 1e-6)
        self.diff2 = Diff("diff2", self.ssys, self.spec_b, 1e-6)
    
    def tearDown(self):
        self.model = None
        self.ssys = None
        self.spec_a = None
        self.spec_b = None
        self.diff1 = None
        self.diff2 = None
    
    def testGetDiff(self):
        self.assertEqual(self.ssys.getDiff("diff1").this, self.diff1.this)
    
    def testDelDiff(self):
        self.ssys.delDiff("diff1")
        self.assertRaises(Exception, self.ssys.getDiff, "diff1")
    
    def testGetAllDiffs(self):
        stored_diffs = [diff.this for diff in self.ssys.getAllDiffs()]
        expected_diffs = [self.diff1.this, self.diff2.this]
        self.assertCountEqual(stored_diffs, expected_diffs)

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(SurfsysCreationTestCase, "test"))
    all_tests.append(unittest.makeSuite(SurfsysPropertyTestCase, "test"))
    all_tests.append(unittest.makeSuite(SurfsysSReacTestCase, "test"))
    all_tests.append(unittest.makeSuite(SurfsysDiffTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
