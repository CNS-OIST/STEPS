########################################################################

# Stochastic production and degradation well-mixed reactions.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Production-Degradation model 
# (see validation/masteq.py)
  
########################################################################

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

print("Reaction - Production and degradation:")
from . import masteq_cp

########################################################################

KCST_f = 100/(Avogadro*1.0e-15) 			# The reaction constant, production
KCST_b = 10			# The reaction constant, degradation
VOL = 1.0e-18

DT = 0.1			# Sampling time-step
INT = 200000.1 		# Sim endtime

# Tolerance for the comparison:
# In tests with good code <1% fail with tolerance of 1.5%
tolerance = 1.5/100  

########################################################################

def test_masteq():
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

    rng = RNG('mt19937', 1000, int(time.time()%4294967295))

    sim = Simulation('Wmdirect', mdl, geom, rng)

    rs = ResultSelector(sim)

    res = rs.comp1.SA.Count

    sim.toSave(res, dt=DT)

    sim.newRun()
    sim.restore('./validation_cp/cp/masteq')
    sim.run(INT)

    def fact(x): return (1 if x==0 else x * fact(x-1))

    # Do cumulative count, but not comparing them all. 
    # Don't get over 50 (I hope)
    steps_n_res = numpy.zeros(50)
    ntpnts = len(res.time[0])
    for r in res.data[0,:,0]: steps_n_res[int(r)]+=1
    for s in range(50): steps_n_res[s] = steps_n_res[s]/ntpnts

    passed = True
    max_err = 0.0

    k1 = KCST_b
    k2 = KCST_f*(Avogadro*1.0e-15)

    # Compare 5 to 15
    for m in range(5, 16):
        analy = (1.0/fact(m))*math.pow((k2/k1), m)*math.exp(-(k2/k1))
        assert(tol_funcs.tolerable(steps_n_res[m], analy, tolerance))
        

