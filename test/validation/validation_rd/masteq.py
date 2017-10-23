########################################################################

# Stochastic production and degradation well-mixed reactions.

# AIMS: to verify STEPS well-mixed stochastic solver 'Wmdirect' 
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

import math
import numpy
import time 

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

from . import tol_funcs

########################################################################

def test_masteq():
    "Reaction - Production and degradation (Wmdirect)"

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

    mdl  = smod.Model()

    A = smod.Spec('A', mdl)

    volsys = smod.Volsys('vsys',mdl)

    # Production
    R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A], kcst = KCST_f)
    R2 = smod.Reac('R2', volsys, lhs = [A], rhs = [], kcst = KCST_b)

    geom = sgeom.Geom()

    comp1 = sgeom.Comp('comp1', geom, VOL)
    comp1.addVolsys('vsys')

    rng = srng.create('r123', 1000)
    rng.initialize(1000)

    sim = ssolv.Wmdirect(mdl, geom, rng)
    sim.reset()

    tpnts = numpy.arange(0.0, INT, DT)
    ntpnts = tpnts.shape[0]

    res = numpy.zeros([ntpnts])

    sim.reset()
    sim.setCompCount('comp1', 'A', 0)

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
        assert tol_funcs.tolerable(steps_n_res[m], analy, tolerance)

########################################################################
# END


