########################################################################

# Stochastic second-order irreversible reaction: [A]0=[B]0
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Reaction model (see validation/second_order_irev_AA.py)
  
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
        rng.initialize(int(time.time()%4294967295))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = numpy.zeros([NITER, ntpnts, 3])

        new_dir = './validation_cp/cp/'
        os.makedirs(new_dir, exist_ok=True)

        sim.reset()
        sim.setCompSpecConc('comp1', 'A', CONCA)
        sim.setCompSpecConc('comp1', 'B', CONCB)
        sim.checkpoint('./validation_cp/cp/second_order_irev_AA')


    def test_soirevAA(self):
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
        rng.initialize(int(time.time()%4294967295))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = numpy.zeros([NITER, ntpnts, 3])

        seed = int(time.time()%4294967295)
        for i in range (0, NITER):
            sim.restore('./validation_cp/cp/second_order_irev_AA')
            rng.initialize(seed)
            seed += 1
            for t in range(0, ntpnts):
                sim.run(tpnts[t])
                res_m[i, t, 0] = sim.getCompSpecConc('comp1', 'A')
                res_m[i, t, 1] = sim.getCompSpecConc('comp1', 'B')

        mean_res = numpy.mean(res_m, 0)

        invA = numpy.zeros(ntpnts)
        invB = numpy.zeros(ntpnts)
        lineA  = numpy.zeros(ntpnts)
        lineB = numpy.zeros(ntpnts)

        max_err=0.0
        for i in range(ntpnts):
            invA[i] = (1.0/mean_res[i][0])
            invB[i] = (1.0/mean_res[i][1])
            lineA[i] = (1.0/CONCA +((tpnts[i]*KCST)))
            lineB[i] = (1.0/CONCB + ((tpnts[i]*KCST)))
            
            self.assertTrue(tol_funcs.tolerable(invA[i], lineA[i], tolerance))
            self.assertTrue(tol_funcs.tolerable(invB[i], lineB[i], tolerance))
            

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestSecondOrderIrevAA))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

