# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Sets system path and provides common testing tools.
# Call before other (STEPS-specific) import statements.
import tools
tools.add_steps_path()

import unittest
import steps.model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class ModelTestCase(unittest.TestCase):

    """Collection of unit tests for basic operations on a kinetic model.
    """

    def setUp(self):
        self.m = steps.model.Model()
        self.x1 = steps.model.Species('x1', self.m)
        self.x2 = steps.model.Species('x2', self.m)
        self.y = steps.model.Species('y', self.m)
        self.v1 = steps.model.Volsys('v1', self.m)
        self.v2 = steps.model.Volsys('v2', self.m)

    def tearDown(self):
        pass

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def testCreateDuplicateSpecies(self):
        # Creating species with an id possessed by other species must fail.
        try: x1b = steps.model.Species('x1', self.m)
        except: pass
        else: self.fail('Creating species with duplicate id should fail')
        # Creating species with id possessed by volsys shouldn't fail.
        try: x1c = steps.model.Species('v1', self.m)
        except: self.fail('Creating species sharing name with volsys failed')

    def testCreateDuplicateVolsys(self):
        # Creating volsys with id possessed by volsys shouldn't fail.
        try: x1 = steps.model.Volsys('x1', self.m)
        except: self.fail('Creation volsys sharing name with species failed')
        # Creating volsys with an id possessed by other volsys must fail.
        try: v1b = steps.model.Volsys('v1', self.m)
        except: pass
        else: self.fail('Creating volsys with duplicate id should fail')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == '__main__': unittest.main()

# END
