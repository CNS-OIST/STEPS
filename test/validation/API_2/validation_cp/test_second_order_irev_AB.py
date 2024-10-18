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

# Stochastic second-order irreversible reaction: [A]0 != [B]0
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Reaction model (see validation/second_order_irev_AB.py)
  
########################################################################

import unittest

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import math
import numpy
import os
import random
import time 

from . import tol_funcs

########################################################################

KCST = 5.0e6			# The reaction constant

CONCA = 1.0e-6
n = 2
CONCB = CONCA/n
VOL = 9.0e-18			

NITER = 1000			# The number of iterations
DT = 0.1			# Sampling time-step
INT = 1.0			# Sim endtime

# In test runs, with good code, <0.1% will fail with a tolerance of 1% 
tolerance = 1.0/100

########################################################################

class TestSecondOrderIrevAB(unittest.TestCase):

    def setUp(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA, SB, SC = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                SA + SB >r['R1']> SC
                r['R1'].K = KCST

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('mt19937', 512, int(random.random()*1000))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        sim.newRun()
        sim.comp1.SA.Conc = CONCA
        sim.comp1.SB.Conc = CONCB

        new_dir = './validation_cp/cp/'
        os.makedirs(new_dir, exist_ok=True)

        sim.checkpoint('./validation_cp/cp/second_order_irev_AB')

    def test_soirevAB(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA, SB, SC = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                SA + SB >r['R1']> SC
                r['R1'].K = KCST

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('mt19937', 512, int(random.random()*1000))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        rs = ResultSelector(sim)

        res = rs.comp1.LIST(SA, SB).Conc

        sim.toSave(res, dt=DT)

        seed = time.time()%4294967295
        for i in range (0, NITER):
            sim.newRun()
            sim.restore('./validation_cp/cp/second_order_irev_AB')
            rng.initialize(seed)
            seed += 1
            sim.run(INT)

        mean_res = numpy.mean(res.data, 0)

        SC = CONCA-CONCB
        for t, (SA, SB) in zip(res.time[0], mean_res):
            lnBA = math.log(SB/SA)
            lineAB = math.log(CONCB/CONCA) -SC*KCST*t
            
            self.assertTrue(tol_funcs.tolerable(lnBA, lineAB, tolerance))


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestSecondOrderIrevAB))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
