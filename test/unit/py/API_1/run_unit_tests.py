# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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

import sys

import nose
import unittest

import model_test
import directional_dcst_test
import sdiff_bugfix_test
import getROIArea_bugfix_test
import setCompVol_bugfix_test
import multi_sys_test
import gmsh_multitag_test
import genPointsInTets_bugfix_test
import missing_solver_methods_test
import utilities_meshctrl_test
import tetODE_setPatchSReacK_bugfix_test

def suite():
    all_tests = [
        model_test.suite(), directional_dcst_test.suite(), sdiff_bugfix_test.suite(),
        getROIArea_bugfix_test.suite(), multi_sys_test.suite(),
        setCompVol_bugfix_test.suite(), gmsh_multitag_test.suite(),
        genPointsInTets_bugfix_test.suite(), missing_solver_methods_test.suite(),
        utilities_meshctrl_test.suite(),
        tetODE_setPatchSReacK_bugfix_test.suite(),
    ]
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2, buffer=True).run(suite())
    if not result.wasSuccessful():
        sys.exit(1)
    nose.run(argv=['test.py', '-v', '-w', 'distribution_test', '--all-modules'])
