####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

########################################################################

# Stochastic second-order irreversible reaction: [A]0=[B]0
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Reaction model (see validation/second_order_irev_AA.py)
  
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
import time 

from . import tol_funcs

########################################################################
KCST = 50.0e6			# The reaction constant

CONCA = 20.0e-6
CONCB = CONCA

VOL = 9.0e-18			

NITER = 1000			# The number of iterations
DT = 0.1                # Sampling time-step
INT = 1.1               # Sim endtime

# In test runs, with good code, <0.1% will fail with a tolerance of 1% 
tolerance = 1.0/100

########################################################################

class TestSecondOrderIrevAA(unittest.TestCase):

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

        rng = RNG('mt19937', 512, int(time.time()%4294967295))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        sim.newRun()
        sim.comp1.SA.Conc = CONCA
        sim.comp1.SB.Conc = CONCB

        new_dir = './validation_cp/cp/'
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)

        sim.checkpoint('./validation_cp/cp/second_order_irev_AA')

    def test_soirevAA(self):
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

        rng = RNG('mt19937', 512, int(time.time()%4294967295))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        rs = ResultSelector(sim)

        res = rs.comp1.LIST(SA, SB).Conc

        sim.toSave(res, dt=DT)

        for i in range (0, NITER):
            sim.newRun()
            sim.restore('./validation_cp/cp/second_order_irev_AA')
            sim.run(INT)

        mean_res = numpy.mean(res.data, 0)

        for t, (CA, CB) in zip(res.time[0], mean_res):
            invA = (1.0/CA)
            invB = (1.0/CB)
            lineA = (1.0/CONCA +((t*KCST)))
            lineB = (1.0/CONCB + ((t*KCST)))
            
            self.assertTrue(tol_funcs.tolerable(invA, lineA, tolerance))
            self.assertTrue(tol_funcs.tolerable(invB, lineB, tolerance))
            

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestSecondOrderIrevAA, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
