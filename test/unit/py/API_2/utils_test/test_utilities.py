####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

""" Unit tests for utilities in utils.py."""

import unittest

from steps import interface

import steps.utils as utils
import steps.model as model

class VariousUtilities(unittest.TestCase):
    """Test units declaration and conversion."""

    def testFreezeAfterInit(self):
        class Control:
            def __init__(self):
                self.a = 5

        @utils.FreezeAfterInit
        class Freeze:
            def __init__(self):
                self.a = 5

        @utils.FreezeAfterInit
        class FreezeSubclass(Freeze):
            def __init__(self):
                super().__init__()
                self.b = 6

        c = Control()
        self.assertEqual(c.a, 5)
        c.A = 10
        self.assertTrue(hasattr(c, 'A'))
        self.assertEqual(c.A, 10)

        f = Freeze()
        self.assertEqual(f.a, 5)
        with self.assertRaises(AttributeError):
            f.A = 10
        self.assertFalse(hasattr(f, 'A'))

        fs = FreezeSubclass()
        self.assertEqual(fs.a, 5)
        self.assertEqual(fs.b, 6)
        with self.assertRaises(AttributeError):
            fs.A = 10
        with self.assertRaises(AttributeError):
            fs.B = 10
        self.assertFalse(hasattr(fs, 'A'))
        self.assertFalse(hasattr(fs, 'B'))

    def testUnusedKwargs(self):
        with self.assertRaises(Exception):
            mdl = model.Model(test=5)
        with self.assertRaises(Exception):
            mdl = model.Model('unused')


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(VariousUtilities, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


