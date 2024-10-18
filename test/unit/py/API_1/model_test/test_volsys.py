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

""" Unit tests for steps.model.Volsys class."""

import unittest
import random
from steps.model import *

class VolsysCreationTestCase(unittest.TestCase):
    """ Tests for volsys construction. """
    def setUp(self):
        self.model = Model()
    
    def tearDown(self):
        self.model = None
    
    def testCreateVolsys(self):
        vsys = Volsys("vsys", self.model)
        self.assertIsInstance(vsys, Volsys)


class VolsysPropertyTestCase(unittest.TestCase):
    """ Property tests for volsys. """
    def setUp(self):
        self.model = Model()
        self.vsys = Volsys("vsys", self.model)
    
    def tearDown(self):
        self.model = None
        self.vsys = None
    
    def testSetGetID(self):
        self.assertEqual(self.vsys.getID(), "vsys")
        self.vsys.setID("new_id")
        self.assertEqual(self.vsys.getID(), "new_id")

    def testGetModel(self):
        self.assertEqual(self.vsys.getModel().this, self.model.this)

class VolsysReacTestCase(unittest.TestCase):
    """ Reac tests for volsys. """
    def setUp(self):
        self.model = Model()
        self.vsys = Volsys("vsys", self.model)
        self.spec_a = Spec("A", self.model)
        self.spec_b = Spec("B", self.model)
        self.spec_c = Spec("C", self.model)
        self.reac1 = Reac("reac1", self.vsys, [self.spec_a, self.spec_b], [self.spec_c], 1e-6)
        self.reac2 = Reac("reac2", self.vsys, [self.spec_a, self.spec_b], [self.spec_c], 1e-6)
    
    def tearDown(self):
        self.model = None
        self.vsys = None
        self.spec_a = None
        self.spec_b = None
        self.spec_c = None
        self.reac1 = None
        self.reac2 = None
    
    def testGetReac(self):
        self.assertEqual(self.vsys.getReac("reac1").this, self.reac1.this)
    
    def testDelReac(self):
        self.vsys.delReac("reac1")
        self.assertRaises(Exception, self.vsys.getReac, "reac1")

    def testGetAllReacs(self):
        stored_reacs = [reac.this for reac in self.vsys.getAllReacs()]
        expected_reacs = [self.reac1.this, self.reac2.this]
        self.assertCountEqual(stored_reacs, expected_reacs)

class VolsysDiffTestCase(unittest.TestCase):
    """ Diff tests for volsys. """
    def setUp(self):
        self.model = Model()
        self.vsys = Volsys("vsys", self.model)
        self.spec_a = Spec("A", self.model)
        self.spec_b = Spec("B", self.model)
        self.diff1 = Diff("diff1", self.vsys, self.spec_a, 1e-6)
        self.diff2 = Diff("diff2", self.vsys, self.spec_b, 1e-6)
    
    def tearDown(self):
        self.model = None
        self.vsys = None
        self.spec_a = None
        self.spec_b = None
        self.diff1 = None
        self.diff2 = None
    
    def testGetDiff(self):
        self.assertEqual(self.vsys.getDiff("diff1").this, self.diff1.this)
    
    def testDelDiff(self):
        self.vsys.delDiff("diff1")
        self.assertRaises(Exception, self.vsys.getDiff, "diff1")
    
    def testGetAllDiffs(self):
        stored_diffs = [diff.this for diff in self.vsys.getAllDiffs()]
        expected_diffs = [self.diff1.this, self.diff2.this]
        self.assertCountEqual(stored_diffs, expected_diffs)

class VolsysSpecTestCase(unittest.TestCase):
    """ Spec tests for volsys. """

    def setUp(self):
        self.model = Model()
        self.vsys = Volsys("vsys", self.model)
        self.spec_a = Spec("A", self.model)
        self.spec_b = Spec("B", self.model)
        self.spec_c = Spec("C", self.model)
        self.diff1 = Diff("diff1", self.vsys, self.spec_a, 1e-6)
        self.diff2 = Diff("diff2", self.vsys, self.spec_b, 1e-6)
    
    def tearDown(self):
        self.model = None
        self.vsys = None
        self.spec_a = None
        self.spec_b = None
    
    def testGetAllSpecs(self):
        stored_specs = [spec.this for spec in self.vsys.getAllSpecs()]
        expected_specs = [self.spec_a.this, self.spec_b.this]
        self.assertCountEqual(stored_specs, expected_specs)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(VolsysCreationTestCase))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(VolsysPropertyTestCase))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(VolsysReacTestCase))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(VolsysDiffTestCase))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(VolsysSpecTestCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
