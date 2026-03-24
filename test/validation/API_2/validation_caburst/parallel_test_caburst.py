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
import numpy as np
import pickle

import os.path as path
from . import caBurstFullModel

FILEDIR = path.dirname(path.abspath(__file__))

class TestCaburst(unittest.TestCase):

    def setUp(self):
        global C

        # defaults
        C={ 'meshdir': path.join(FILEDIR, 'meshes', 'split_8'),
            'mesh': 'pyramid_1_0.2_112tets',
            'datadir': path.join(FILEDIR, 'data'),
            'v1_3': 'v1_3.txt', 
            'v2_3': 'v2_3.txt', 
            'v3_3': 'v3_3.txt', 
            'v4_3': 'v4_3.txt',
            'v1_4': 'v1_4.txt', 
            'v2_4': 'v2_4.txt', 
            'v3_4': 'v3_4.txt', 
            'v4_4': 'v4_4.txt' }

    def test_caburst_n8(self):

        seed = 123
        niter = 10

        meshfile = path.join(C['meshdir'], C['mesh'])

        sim3data, sim3time = caBurstFullModel.run(seed, meshfile, 3, niter)
        sim4data, sim4time = caBurstFullModel.run(seed, meshfile, 4, niter)

        sim3data_mean = np.mean(sim3data, axis=0)
        sim4data_mean = np.mean(sim4data, axis=0)

        for ref_version in ['3', '4']:
            v1_ifile = open(path.join(C['datadir'], C[f'v1_{ref_version}']), 'rb')
            v2_ifile = open(path.join(C['datadir'], C[f'v2_{ref_version}']), 'rb')
            v3_ifile = open(path.join(C['datadir'], C[f'v3_{ref_version}']), 'rb')
            v4_ifile = open(path.join(C['datadir'], C[f'v4_{ref_version}']), 'rb')

            v1_mean, v1_std = pickle.load(v1_ifile), pickle.load(v1_ifile)
            v2_mean, v2_std = pickle.load(v2_ifile), pickle.load(v2_ifile)
            v3_mean, v3_std = pickle.load(v3_ifile), pickle.load(v3_ifile)
            v4_mean, v4_std = pickle.load(v4_ifile), pickle.load(v4_ifile)

            v1_ifile.close()
            v2_ifile.close()
            v3_ifile.close()
            v4_ifile.close()

            sigma = 1.5
            lb = 50 # Remove first 10ms, pre-stimulus, which enables a lower sigma
            
            self.assertTrue((sim3data_mean[:, 0][lb:] < (v1_mean+sigma*v1_std)[lb:]).all())
            self.assertTrue((sim3data_mean[:, 1][lb:] < (v2_mean+sigma*v2_std)[lb:]).all())
            self.assertTrue((sim3data_mean[:, 2][lb:] < (v3_mean+sigma*v3_std)[lb:]).all())
            self.assertTrue((sim3data_mean[:, 3][lb:] < (v4_mean+sigma*v4_std)[lb:]).all())

            self.assertTrue((sim3data_mean[:, 0][lb:] > (v1_mean-sigma*v1_std)[lb:]).all())
            self.assertTrue((sim3data_mean[:, 1][lb:] > (v2_mean-sigma*v2_std)[lb:]).all())
            self.assertTrue((sim3data_mean[:, 2][lb:] > (v3_mean-sigma*v3_std)[lb:]).all())

            self.assertTrue((sim4data_mean[:, 0][lb:] < (v1_mean+sigma*v1_std)[lb:]).all())
            self.assertTrue((sim4data_mean[:, 1][lb:] < (v2_mean+sigma*v2_std)[lb:]).all())
            self.assertTrue((sim4data_mean[:, 2][lb:] < (v3_mean+sigma*v3_std)[lb:]).all())
            self.assertTrue((sim4data_mean[:, 3][lb:] < (v4_mean+sigma*v4_std)[lb:]).all())
            
            self.assertTrue((sim4data_mean[:, 0][lb:] > (v1_mean-sigma*v1_std)[lb:]).all())
            self.assertTrue((sim4data_mean[:, 1][lb:] > (v2_mean-sigma*v2_std)[lb:]).all())
            self.assertTrue((sim4data_mean[:, 2][lb:] > (v3_mean-sigma*v3_std)[lb:]).all())

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestCaburst))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
