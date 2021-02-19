########################################################################

# Stochastic first-order reversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Reversible 
# Reaction model (see validation/first_order_rev.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import time 
import numpy

from scipy.constants import Avogadro
from . import tol_funcs

print("Reaction - First order, reversible:")
from . import first_order_rev_cp

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

def test_forev():
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

    for i in range (0, NITER):
        sim.newRun()
        sim.restore('./validation_cp/cp/first_order_rev')
        sim.run(INT)

    mean_res = numpy.mean(res.data, 0) * 1e6

    Aeq = COUNT*(KCST_b/KCST_f)/(1+(KCST_b/KCST_f))/(VOL*Avogadro * 1e3)*1e6
    Beq = (COUNT/(VOL*Avogadro * 1e3))*1e6 -Aeq

    max_err = 0.0
    passed = True
    for i in range(len(res.time[0])):
        if i < 7:
            continue
        assert(tol_funcs.tolerable(mean_res[i,0], Aeq, tolerance))
        assert(tol_funcs.tolerable(mean_res[i,1], Beq, tolerance))


