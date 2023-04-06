########################################################################

# Stochastic production and degradation well-mixed reactions.

# AIMS: to verify STEPS well-mixed stochastic solver 'Wmrssa' 
# solves the master equation correctly when applied to a  
# zero-order production and first-order degradation process.

# For a more detailed description of the analytical system and 
# equivalent STEPS model see:
# http://www.biomedcentral.com/content/supplementary/1752-0509-6-36-s4.pdf
# "Production and degradation reactions"

# Verification also takes place of the necessary steps to build the model, 
# such as well-mixed compartment creation, random-number generator 
# construction and initialization, and recording from a well-mixed 
# compartment. 

# A 1.5% tolerance is imposed when comparing the stationary distribution 
# from 200000s of STEPS stochastic simulation to the analytical solution 
# to the chemical master equation, in the range 5-15 molecules. 
# There is an expected probability of failure of < 1%.
  
########################################################################

import unittest

import steps.interface

import math
import numpy
import time 
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

from scipy.constants import Avogadro
from . import tol_funcs

########################################################################

class TestRDMasteqRSSA(unittest.TestCase):
    def test_masteq_rssa(self):
        "Reaction - Production and degradation (Wmrssa)"

        ########################################################################

        KCST_f = 100/(Avogadro*1.0e-15) 			# The reaction constant, production
        KCST_b = 10			# The reaction constant, degradation
        VOL = 1.0e-18

        DT = 0.1			# Sampling time-step
        INT = 200000.0 		# Sim endtime

        # Tolerance for the comparison:
        # In tests with good code <1% fail with tolerance of 1.5%
        tolerance = 1.5/100  

        ########################################################################
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                # Production
                None <r['R1']> SA
                r['R1'].K = KCST_f, KCST_b

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(volsys, VOL)

        rng = RNG('r123', 1000, 1000)

        sim = Simulation('Wmrssa', mdl, geom, rng)

        rs = ResultSelector(sim)

        res = rs.comp1.SA.Count

        sim.toSave(res, dt=DT)

        sim.newRun()
        sim.comp1.SA.Count = 0
        sim.run(INT)

        def fact(x): return (1 if x==0 else x * fact(x-1))

        # Do cumulative count, but not comparing them all. 
        # Don't get over 50 (I hope)
        steps_n_res = numpy.zeros(50)
        for r in res.data[0,:,0]: steps_n_res[int(r)]+=1
        for s in range(50): steps_n_res[s] = steps_n_res[s]/len(res.time[0])

        passed = True
        max_err = 0.0

        k1 = KCST_b
        k2 = KCST_f*(Avogadro*1.0e-15)

        # Compare 5 to 15
        for m in range(5, 16):
            analy = (1.0/fact(m))*math.pow((k2/k1), m)*math.exp(-(k2/k1))
            self.assertTrue(tol_funcs.tolerable(steps_n_res[m], analy, tolerance))

    ########################################################################

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestRDMasteqRSSA, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
