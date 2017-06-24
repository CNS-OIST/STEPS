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

""" Unit tests for steps.model.Chan class."""

import unittest2
import random
from steps.model import *

class ChanCreationTestCase(unittest2.TestCase):
    """ Tests for channel construction. """
    def setUp(self):
        self.model = Model()
    
    def tearDown(self):
        self.model = None
    
    def testCreateChan(self):
        chan = Chan("A", self.model)
        self.assertIsInstance(chan, Chan)

class ChanPropertyTestCase(unittest2.TestCase):
    """ Property tests for channel. """
    def setUp(self):
        self.model = Model()
        self.chan = Chan("A", self.model)
    
    def tearDown(self):
        self.model = None
        self.Chan_normal = None
        self.Chan_valence = None
    
    def testSetGetChanID(self):
        self.assertEqual(self.chan.getID(), "A")
        self.chan.setID("new_id")
        self.assertEqual(self.chan.getID(), "new_id")

    def testGetModel(self):
        self.assertEqual(self.chan.getModel().this, self.model.this)

class ChanStateTestCase(unittest2.TestCase):
    """ ChanState tests. """
    def setUp(self):
        self.model = Model()
        self.chan = Chan("chan", self.model)
        self.state1 = ChanState("state1", self.model, self.chan)
        self.state2 = ChanState("state2", self.model, self.chan)
    
    def tearDown(self):
        self.model = None
        self.chan = None
        self.state1 = None
        self.state2 = None

    def testGetChanState(self):
        self.assertEqual(self.chan.getChanState("state1").this, self.state1.this)

    def testGetAllChanStates(self):
        stored_states = [state.this for state in self.chan.getAllChanStates()]
        expected_states = [self.state1.this, self.state2.this]
        self.assertItemsEqual(stored_states, expected_states)

def suite():
    all_tests = []
    all_tests.append(unittest2.makeSuite(ChanCreationTestCase, "test"))
    all_tests.append(unittest2.makeSuite(ChanPropertyTestCase, "test"))
    all_tests.append(unittest2.makeSuite(ChanStateTestCase, "test"))
    return unittest2.TestSuite(all_tests)

if __name__ == "__main__":
    unittest2.TextTestRunner(verbosity=2).run(suite())
