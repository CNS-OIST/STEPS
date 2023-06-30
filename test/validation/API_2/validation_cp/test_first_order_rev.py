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

# Stochastic first-order reversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Reversible 
# Reaction model (see validation/first_order_rev.py)
  
########################################################################

import unittest

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import time 
import os
import numpy

from scipy.constants import Avogadro
from . import tol_funcs

########################################################################

KCST_f = 10.0			# The reaction constant
KCST_b = 2.0

COUNT = 100000				# Can set count or conc
VOL = 6.0e-18

NITER = 10			# The number of iterations
DT = 0.1		# Sampling time-step
INT = 1.0			# Sim endtime

# In test runs, with good code, <0.1% will fail with a tolerance of 1% 
tolerance = 1.0/100

########################################################################

class TestFirstOrderRev(unittest.TestCase):

    def setUp(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA, SB = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                SA <r['R1']> SB
                r['R1'].K = KCST_f, KCST_b

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('mt19937', 512, int(time.time()%4294967295))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        sim.newRun()
        sim.comp1.SA.Count = COUNT
        sim.comp1.SB.Count = 0.0

        new_dir = './validation_cp/cp/'
        os.makedirs(new_dir, exist_ok=True)

        sim.checkpoint('./validation_cp/cp/first_order_rev')

    def test_forev(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA, SB = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                SA <r['R1']> SB
                r['R1'].K = KCST_f, KCST_b

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('mt19937', 512, int(time.time()%4294967295))

        sim = Simulation('Wmdirect', mdl, geom, rng)

        rs = ResultSelector(sim)

        res = rs.comp1.LIST(SA, SB).Conc

        sim.toSave(res, dt=DT)

        seed = time.time()%4294967295
        for i in range (0, NITER):
            sim.newRun()
            sim.restore('./validation_cp/cp/first_order_rev')
            rng.initialize(seed)
            seed += 1
            sim.run(INT)

        mean_res = numpy.mean(res.data, 0) * 1e6

        Aeq = COUNT*(KCST_b/KCST_f)/(1+(KCST_b/KCST_f))/(VOL*Avogadro * 1e3)*1e6
        Beq = (COUNT/(VOL*Avogadro * 1e3))*1e6 -Aeq

        max_err = 0.0
        for i in range(len(res.time[0])):
            if i < 7:
                continue
            self.assertTrue(tol_funcs.tolerable(mean_res[i,0], Aeq, tolerance))
            self.assertTrue(tol_funcs.tolerable(mean_res[i,1], Beq, tolerance))

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestFirstOrderRev, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
