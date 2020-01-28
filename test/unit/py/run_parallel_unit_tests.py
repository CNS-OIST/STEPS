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

import unittest2
import nose
import parallel_diff_sel_test
import parallel_setget_count_test
def suite():
    all_tests = [ parallel_diff_sel_test.suite(), parallel_setget_count_test.suite()]
    return unittest2.TestSuite(all_tests)

if __name__ == "__main__":
    unittest2.TextTestRunner(verbosity=2, buffer=True).run(suite())
    #nose.run(argv=['test.py', '-v', '-w', 'parallel_diff_sel_test', '--all-modules'])
