####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
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
from . import rallpack1

class TestRallpack1(unittest.TestCase):

    def setUp(self):
        global C

        # defaults
        C={ 'meshdir': 'validation_efield/meshes',
            'mesh': 'axon_cube_L1000um_D866nm_1978tets',
            'datadir': 'validation_efield/data/rallpack1_correct',
            'v0data': 'v0',
            'v1data': 'vx',
            'seed': 7 }


    def test_rallpack1(self):
        params = rallpack1.sim_parameters

        meshfile = path.join(C['meshdir'],C['mesh'])
        v0data = path.join(C['datadir'],C['v0data'])
        v1data = path.join(C['datadir'],C['v1data'])
        seed = C['seed']

        simdata, rms_err_0um, rms_err_1000um = rallpack1.run_comparison(seed, meshfile, v0data, v1data)

        print("rms error at 0um = " + str(rms_err_0um))
        print("rms error at 1000um = " + str(rms_err_1000um))
        
        # import matplotlib.pyplot as plt
        # plt.subplot(211)
        # plt.plot(simdata[0,:], simdata[2,:], 'k-' ,label = 'Correct, 0um', linewidth=3)
        # plt.plot(simdata[0,:], simdata[1,:], 'r--', label = 'STEPS, 0um', linewidth=3)
        # plt.legend(loc='best')
        # plt.ylabel('Potential (mV)')
        # plt.subplot(212)
        # plt.plot(simdata[0,:], simdata[4,:], 'k-' ,label = 'Correct, 1000um', linewidth=3)
        # plt.plot(simdata[0,:], simdata[3,:], 'r--', label = 'STEPS, 1000um', linewidth=3)
        # plt.legend(loc='best')
        # plt.ylabel('Potential (mV)')
        # plt.xlabel('Time (ms)')
        # plt.show(block=True)

        max_rms_err = 1.e-3
        self.assertTrue(rms_err_0um < max_rms_err)
        self.assertTrue(rms_err_1000um < max_rms_err)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestRallpack1))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
