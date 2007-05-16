# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Sets system path and provides common testing tools.
# Call before other (STEPS-specific) import statements.
import tools
tools.add_steps_path()

import steps
import unittest

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class IDTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

    def testValidIDs(self):
        """Test the range of valid component id's."""
        self.assertEqual(steps.isValidID('abaca'), True, \
                '\'abaca\' should be a valid id')
        self.assertEqual(steps.isValidID('abaca00'), True, \
                '\'abaca00\' should be a valid id')
        self.assertEqual(steps.isValidID('_a'), True, \
                '\'_a\' should be a valid id')
        self.assertEqual(steps.isValidID('a_00d'), True, \
                '\'a_00d\' should be a valid id')
        self.assertEqual(steps.isValidID(''), False, \
                'An empty string should not be a valid id')
        self.assertEqual(steps.isValidID('00d'), False, \
                '\'00d\' should not be a valid id')
        self.assertEqual(steps.isValidID('  d'), False, \
                '\'  d\' should not be a valid id')
        self.assertEqual(steps.isValidID('d  '), False, \
                '\'d  \' should not be a valid id')
        self.assertRaises(TypeError, steps.isValidID, 0.0)
        self.assertRaises(TypeError, steps.isValidID, 0)
        self.assertRaises(TypeError, steps.isValidID, ['a'])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == '__main__': unittest.main()

# END
