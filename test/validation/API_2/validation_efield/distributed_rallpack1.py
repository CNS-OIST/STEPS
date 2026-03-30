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

class DistributedRallpack1(TestRallpack1Base):

    common_params = dict(
        solver='DistTetOpSplit',
        mesh='axon_cube_L1000um_D866nm_1135tets.msh',
        sim_end=0.05,
    )

    def test_rallpack1_distTetOpSplit(self):
        self._test_rallpack1(**self.common_params, currType='Ohmic', plot=False)

    def test_rallpack1_vdsr_distTetOpSplit(self):
        self._test_rallpack1(**self.common_params, currType='SReacCharge', plot=False)

    def test_rallpack1_distTetOpSplit_multicond_seq(self):
        self._test_rallpack1(**self.common_params, currType='Ohmic', comp_mode='multi_seq', checks=['end'], plot=False)

    def test_rallpack1_distTetOpSplit_multicond_par(self):
        self._test_rallpack1(**self.common_params, currType='Ohmic', comp_mode='multi_par', plot=False)

    def test_rallpack1_distTetOpSplit_multicond_rand(self):
        self._test_rallpack1(**self.common_params, currType='Ohmic', comp_mode='multi_rand', max_rms_err=1e-2, plot=False)

    # This test is not run in the CI, only there to check that DistTetOpSplit behaves as expected when the conduction
    # volume is split in many compartments.
    # def test_rallpack1_distTetOpSplit_multicond_maxcomps(self):
    #     self._test_rallpack1(**self.common_params, currType='Ohmic', comp_mode='multi_maxcomps', plot=True)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistributedRallpack1))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
