########################################################################

# Stochastic first-order reversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Reversible 
# Reaction model (see validation/first_order_rev.py)
  
########################################################################

import unittest

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import os
import time 
import numpy

from . import tol_funcs

########################################################################

KCST_f = 10.0			# The reaction constant
KCST_b = 2.0

COUNT = 100000				# Can set count or conc
VOL = 6.0e-18

NITER = 10			# The number of iterations
DT = 0.1		# Sampling time-step
INT = 1.1			# Sim endtime

# In test runs, with good code, <0.1% will fail with a tolerance of 1% 
tolerance = 1.0/100

########################################################################

class TestFirstOrderRev(unittest.TestCase):

    def setUp(self):
        mdl  = smod.Model()

        A = smod.Spec('A', mdl)
        B = smod.Spec('B', mdl)

        volsys = smod.Volsys('vsys',mdl)

        R1 = smod.Reac('R1', volsys, lhs = [A], rhs = [B], kcst = KCST_f)
        R2 = smod.Reac('R2', volsys, lhs = [B], rhs = [A], kcst = KCST_b)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom, VOL)
        comp1.addVolsys('vsys')

        rng = srng.create('mt19937', 512)
        rng.initialize(int(time.time()%4294967295))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = numpy.zeros([NITER, ntpnts, 2]) 

        new_dir = './validation_cp/cp/'
        os.makedirs(new_dir, exist_ok=True)

        sim.reset()
        sim.setCompSpecCount('comp1', 'A', COUNT)
        sim.setCompSpecCount('comp1', 'B', 0.0)
        sim.checkpoint('./validation_cp/cp/first_order_rev')

    def test_forev(self):
        mdl  = smod.Model()

        A = smod.Spec('A', mdl)
        B = smod.Spec('B', mdl)

        volsys = smod.Volsys('vsys',mdl)

        R1 = smod.Reac('R1', volsys, lhs = [A], rhs = [B], kcst = KCST_f)
        R2 = smod.Reac('R2', volsys, lhs = [B], rhs = [A], kcst = KCST_b)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom, VOL)
        comp1.addVolsys('vsys')

        rng = srng.create('mt19937', 512)
        rng.initialize(int(time.time()%4294967295))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res_m = numpy.zeros([NITER, ntpnts, 2]) 

        seed = int(time.time()%4294967295)
        for i in range (0, NITER):
            sim.restore('./validation_cp/cp/first_order_rev')
            rng.initialize(seed)
            seed += 1
            for t in range(0, ntpnts):
                sim.run(tpnts[t])
                res_m[i, t, 0] = sim.getCompSpecConc('comp1', 'A')*1e6
                res_m[i, t, 1] = sim.getCompSpecConc('comp1', 'B')*1e6

        mean_res = numpy.mean(res_m, 0)

        Aeq = COUNT*(KCST_b/KCST_f)/(1+(KCST_b/KCST_f))/(VOL*6.0221415e26)*1e6
        Beq = (COUNT/(VOL*6.0221415e26))*1e6 -Aeq

        max_err = 0.0
        for i in range(ntpnts):
            if i < 7: continue
            self.assertTrue(tol_funcs.tolerable(mean_res[i,0], Aeq, tolerance))
            self.assertTrue(tol_funcs.tolerable(mean_res[i,1], Beq, tolerance))


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestFirstOrderRev, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

