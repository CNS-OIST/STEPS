########################################################################

# Stochastic second-order irreversible reaction: [A]0 != [B]0
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Reaction model (see validation/second_order_irev_AB.py)
  
########################################################################

import unittest

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import os
import math
import time 
import numpy

from . import tol_funcs

########################################################################

KCST = 5.0e6			# The reaction constant

CONCA = 1.0e-6
n = 2
CONCB = CONCA/n
VOL = 9.0e-18			

NITER = 1000			# The number of iterations
DT = 0.1			# Sampling time-step
INT = 1.1			# Sim endtime

# In test runs, with good code, <0.1% will fail with a tolerance of 1% 
tolerance = 1.0/100

########################################################################

class TestSecondOrderIrevAB(unittest.TestCase):

    def setUp(self):
        mdl  = smod.Model()

        A = smod.Spec('A', mdl)
        B = smod.Spec('B', mdl)
        C = smod.Spec('C', mdl)

        volsys = smod.Volsys('vsys',mdl)

        R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [C], kcst = KCST)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom, VOL)
        comp1.addVolsys('vsys')

        rng = srng.create('mt19937', 512)
        import random
        rng.initialize(int(random.random()*1000))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = numpy.zeros([NITER, ntpnts, 3])

        new_dir = './validation_cp/cp/'
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)

        sim.reset()
        sim.setCompConc('comp1', 'A', CONCA)
        sim.setCompConc('comp1', 'B', CONCB)
        sim.checkpoint('./validation_cp/cp/second_order_irev_AB')


    def test_soirevAB(self):
        mdl  = smod.Model()

        A = smod.Spec('A', mdl)
        B = smod.Spec('B', mdl)
        C = smod.Spec('C', mdl)

        volsys = smod.Volsys('vsys',mdl)

        R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [C], kcst = KCST)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom, VOL)
        comp1.addVolsys('vsys')

        rng = srng.create('mt19937', 512)
        import random
        rng.initialize(int(random.random()*1000))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = numpy.zeros([NITER, ntpnts, 3])

        for i in range (0, NITER):
            sim.restore('./validation_cp/cp/second_order_irev_AB')
            
            for t in range(0, ntpnts):
                sim.run(tpnts[t])
                res_m[i, t, 0] = sim.getCompConc('comp1', 'A')
                res_m[i, t, 1] = sim.getCompConc('comp1', 'B')

        mean_res = numpy.mean(res_m, 0)

        lnBA = numpy.zeros(ntpnts)
        lineAB = numpy.zeros(ntpnts)
        C = CONCA-CONCB
        passed  =True
        max_err = 0.0
        for i in range(ntpnts):
            A = mean_res[i][0]
            B = mean_res[i][1]
            lnBA[i] = math.log(B/A)
            lineAB[i] = math.log(CONCB/CONCA) -C*KCST*tpnts[i]
            
            self.assertTrue(tol_funcs.tolerable(lnBA[i], lineAB[i], tolerance))

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestSecondOrderIrevAB, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

