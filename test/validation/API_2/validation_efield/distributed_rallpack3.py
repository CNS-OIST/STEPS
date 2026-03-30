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

import os.path as path
import unittest

from .test_rallpack3 import TestRallpack3Base


class DistributedRallpack3(TestRallpack3Base):

    common_params = dict(
        mesh="axon_ecs_1000um_1um_1314tets.msh",
        sim_end=0.05,
    )

    def test_rallpack3_distTetOpSplit_ohmic_n2(self):
        self._test_rallpack3(
            **self.common_params, solver="DistTetOpSplit", currType="Ohmic", plot=False
        )

    def test_rallpack3_distTetOpSplit_ghk_n8(self):
        self._test_rallpack3(
            **self.common_params,
            solver="DistTetOpSplit",
            currType="GHK",
            ion_diff=True,
            plot=False,
            SAVE_DT=5e-6,
            v0data=self.V0_GHK_REF,
            v1data=self.V1_GHK_REF,
            max_rms_err=5,
        )

    def test_rallpack3_distTetOpSplit_ghk_rleaping_n8(self):
        self._test_rallpack3(
            **self.common_params,
            solver="RLeapDistTetOpSplit",
            currType="GHK",
            ion_diff=True,
            plot=False,
            SAVE_DT=5e-6,
            v0data=self.V0_GHK_REF,
            v1data=self.V1_GHK_REF,
            max_rms_err=5,
        )


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistributedRallpack3))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
