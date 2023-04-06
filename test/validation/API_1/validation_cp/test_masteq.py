########################################################################

# Stochastic production and degradation well-mixed reactions.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Production-Degradation model 
# (see validation/masteq.py)
  
########################################################################

import unittest

import os
import math
import numpy
import time 

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

from . import tol_funcs

########################################################################

KCST_f = 100/(6.022e23*1.0e-15) 			# The reaction constant, production
KCST_b = 10			# The reaction constant, degradation
VOL = 1.0e-18

DT = 0.1			# Sampling time-step
INT = 200000.1 		# Sim endtime

# Tolerance for the comparison:
# In tests with good code <1% fail with tolerance of 1.5%
tolerance = 1.5/100  

########################################################################

class TestMastEq(unittest.TestCase):

    def setUp(self):

        mdl  = smod.Model()

        A = smod.Spec('A', mdl)

        volsys = smod.Volsys('vsys',mdl)

        # Production
        R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A], kcst = KCST_f)
        R2 = smod.Reac('R2', volsys, lhs = [A], rhs = [], kcst = KCST_b)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom, VOL)
        comp1.addVolsys('vsys')

        rng = srng.create('mt19937', 1000)
        rng.initialize(int(time.time()%4294967295))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res = numpy.zeros([ntpnts])

        sim.reset()
        sim.setCompCount('comp1', 'A', 0)

        new_dir = './validation_cp/cp/'
        os.makedirs(new_dir, exist_ok=True)

        sim.checkpoint('./validation_cp/cp/masteq')


    def test_masteq(self):
        mdl  = smod.Model()

        A = smod.Spec('A', mdl)

        volsys = smod.Volsys('vsys',mdl)

        # Production
        R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A], kcst = KCST_f)
        R2 = smod.Reac('R2', volsys, lhs = [A], rhs = [], kcst = KCST_b)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom, VOL)
        comp1.addVolsys('vsys')

        rng = srng.create('mt19937', 1000)
        rng.initialize(int(time.time()%4294967295))

        sim = ssolv.Wmdirect(mdl, geom, rng)
        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        res = numpy.zeros([ntpnts])

        sim.restore('./validation_cp/cp/masteq')

        for t in range(0, ntpnts):
            sim.run(tpnts[t])
            res[t] = sim.getCompCount('comp1', 'A')

        def fact(x): return (1 if x==0 else x * fact(x-1))

        # Do cumulative count, but not comparing them all. 
        # Don't get over 50 (I hope)
        steps_n_res = numpy.zeros(50)
        for r in res: steps_n_res[int(r)]+=1
        for s in range(50): steps_n_res[s] = steps_n_res[s]/ntpnts

        passed = True
        max_err = 0.0

        k1 = KCST_b
        k2 = KCST_f*(6.022e23*1.0e-15)

        # Compare 5 to 15
        for m in range(5, 16):
            analy = (1.0/fact(m))*math.pow((k2/k1), m)*math.exp(-(k2/k1))
            self.assertTrue(tol_funcs.tolerable(steps_n_res[m], analy, tolerance))
            

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestMastEq, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

