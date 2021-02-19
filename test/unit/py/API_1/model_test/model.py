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

""" Unit tests for steps.model.Model class."""

import unittest
import random
from steps.model import *

class ModelCreationTestCase(unittest.TestCase):
    """ Tests for model construction. """
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def testCreateModel(self):
        model = Model()
        self.assertIsInstance(model, Model)

class ModelSpecTestCase(unittest.TestCase):
    """   Tests for species data access. """
    
    def setUp(self):
        self.model = Model()
        self.spec_a = Spec("A", self.model)
        self.spec_b = Spec("B", self.model)
    
    def tearDown(self):
        self.model = None
        self.spec_a = None
        self.spec_b = None
    
    def testGetSpec(self):
        self.assertEqual(self.model.getSpec("A").this, self.spec_a.this)

    def testDelSpec(self):
        self.model.delSpec("B")
        self.assertRaises(Exception, self.model.getSpec, "B")

    def testGetAllSpecs(self):
        stored_specs = [spec.this for spec in self.model.getAllSpecs()]
        expected_specs = [self.spec_a.this, self.spec_b.this]
        
        self.assertCountEqual(stored_specs, expected_specs)

class ModelChanTestCase(unittest.TestCase):
    """   Tests for channel data access. """
    
    def setUp(self):
        self.model = Model()
        self.chan_a = Chan("A", self.model)
        self.chan_b = Chan("B", self.model)
    
    def tearDown(self):
        self.model = None
        self.chan_a = None
        self.chan_b = None
    
    def testGetChan(self):
        self.model.getChan("A")
        self.assertEqual(self.model.getChan("A").this, self.chan_a.this)
    
    def testGetAllChans(self):
        stored_chans = [chan.this for chan in self.model.getAllChans()]
        expected_chans = [self.chan_a.this, self.chan_b.this]
        
        self.assertCountEqual(stored_chans, expected_chans)
        
class ModelVolsysTestCase(unittest.TestCase):
    """   Tests for volume system data access. """
    
    def setUp(self):
        self.model = Model()
        self.vsys_a = Volsys("A", self.model)
        self.vsys_b = Volsys("B", self.model)
    
    def tearDown(self):
        self.model = None
        self.vsys_a = None
        self.vsys_b = None
    
    def testGetVolsys(self):
        self.assertEqual(self.model.getVolsys("A").this, self.vsys_a.this)
    
    def testDelVolsys(self):
        self.model.delVolsys("B")
        self.assertRaises(Exception, self.model.getVolsys, "B")
    
    def testGetAllVolsyss(self):
        stored_vsyss = [vsys.this for vsys in self.model.getAllVolsyss()]
        expected_vsyss = [self.vsys_a.this, self.vsys_b.this]
        
        self.assertCountEqual(stored_vsyss, expected_vsyss)

class ModelSurfsysTestCase(unittest.TestCase):
    """   Tests for surface system data access. """
    
    def setUp(self):
        self.model = Model()
        self.ssys_a = Surfsys("A", self.model)
        self.ssys_b = Surfsys("B", self.model)
    
    def tearDown(self):
        self.model = None
        self.ssys_a = None
        self.ssys_b = None
    
    def testGetSurfsys(self):
        self.assertEqual(self.model.getSurfsys("A").this, self.ssys_a.this)
    
    def testDelSurfsys(self):
        self.model.delSurfsys("B")
        self.assertRaises(Exception, self.model.getSurfsys, "B")
    
    def testGetAllSurfsyss(self):
        stored_ssyss = [ssys.this for ssys in self.model.getAllSurfsyss()]
        expected_ssyss = [self.ssys_a.this, self.ssys_b.this]
        
        self.assertCountEqual(stored_ssyss, expected_ssyss)

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ModelCreationTestCase, "test"))
    all_tests.append(unittest.makeSuite(ModelSpecTestCase, "test"))
    all_tests.append(unittest.makeSuite(ModelChanTestCase, "test"))
    all_tests.append(unittest.makeSuite(ModelVolsysTestCase, "test"))
    all_tests.append(unittest.makeSuite(ModelSurfsysTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
