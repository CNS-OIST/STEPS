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

########################################################################

# Stochastic first-order irreversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Irreversible 
# Reaction model (see validation/first_order_irev.py)
  
########################################################################

import unittest

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import numpy as np
import os
import time

from . import tol_funcs

########################################################################

KCST = 5			# The reaction constant
N = 50				# Can set count or conc
VOL = 1.0e-18

NITER = 100000		# The number of iterations
DT = 0.1			# Sampling time-step
INT = 1.0			# Sim endtime

# Tolerance for the comparison:
# In test runs, with good code, < 1%  will fail with a 1.5% tolerance
tolerance = 2.0/100  

########################################################################

class TestFirstOrderIrev(unittest.TestCase):

    def setUp(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                SA >r['R1']> None
                r['R1'].K = KCST

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('mt19937', 1000, int(time.time()%4294967295))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        tpnts = np.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = np.zeros([NITER, ntpnts, 1])
        res_std1 = np.zeros([ntpnts, 1])
        res_std2 = np.zeros([ntpnts, 1])

        sim.newRun()
        sim.comp1.SA.Count = N

        new_dir = './validation_cp/cp/'
        os.makedirs(new_dir, exist_ok=True)

        sim.checkpoint('./validation_cp/cp/first_order_irev')

    def test_foirev(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                SA >r['R1']> None
                r['R1'].K = KCST

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('mt19937', 1000, int(time.time()%4294967295))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        rs = ResultSelector(sim)

        res = rs.comp1.SA.Count

        sim.toSave(res, dt=DT)

        seed = time.time()%4294967295
        for i in range (0, NITER):
            sim.newRun()
            sim.restore('./validation_cp/cp/first_order_irev')
            rng.initialize(seed)
            seed += 1
            sim.run(INT)

        mean_res = np.mean(res.data, 0)
        std_res = np.std(res.data, 0)

        m_tol = 0
        s_tol=0

        for i in range(len(res.time[0])):
            if i == 0:
                continue
            analy = N*np.exp(-KCST*res.time[0,i])
            std = np.power((N*(np.exp(-KCST*res.time[0,i]))*(1-(np.exp(-KCST*res.time[0,i])))), 0.5)
            self.assertTrue(tol_funcs.tolerable(analy, mean_res[i], tolerance))
            self.assertTrue(tol_funcs.tolerable(std, std_res[i], tolerance))

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestFirstOrderIrev))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

