####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

import unittest

import os.path as path
from .test_rallpack1 import TestRallpack1Base

class ParallelRallpack1(TestRallpack1Base):

    def test_rallpack1_tetOpSplit(self):
        self._test_rallpack1(solver='TetOpSplit', currType='Ohmic', sim_end=0.05, plot=False)

    def test_rallpack1_vdsr_tetOpSplit(self):
        self._test_rallpack1(solver='TetOpSplit', currType='SReacCharge', sim_end=0.05, plot=False)

    def test_rallpack1_tetVesicle(self):
        self._test_rallpack1(solver='TetVesicle', currType='Ohmic', sim_end=0.05, plot=False)

    def test_rallpack1_vdsr_tetVesicle(self):
        self._test_rallpack1(solver='TetVesicle', currType='SReacCharge', sim_end=0.05, plot=False)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(ParallelRallpack1))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
