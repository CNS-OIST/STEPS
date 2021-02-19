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

import sys
import unittest

import model_test
import complex_test
import geom_test
import sim_test

def suite():
    all_tests = [
        model_test.suite(), complex_test.suite(), geom_test.suite(), sim_test.suite()
    ]
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2, buffer=True).run(suite())
    if not result.wasSuccessful():
        sys.exit(1)

