########################################################################

# Stochastic first-order reversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Reversible 
# Reaction model (see validation/first_order_rev.py)
  
########################################################################

from __future__ import print_function, absolute_import
import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import time 
import numpy

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

    for i in range (0, NITER):
        sim.restore('./validation_cp/cp/first_order_rev')
        for t in range(0, ntpnts):
            sim.run(tpnts[t])
            res_m[i, t, 0] = sim.getCompConc('comp1', 'A')*1e6
            res_m[i, t, 1] = sim.getCompConc('comp1', 'B')*1e6

    mean_res = numpy.mean(res_m, 0)

    Aeq = COUNT*(KCST_b/KCST_f)/(1+(KCST_b/KCST_f))/(VOL*6.0221415e26)*1e6
    Beq = (COUNT/(VOL*6.0221415e26))*1e6 -Aeq

    max_err = 0.0
    passed = True
    for i in range(ntpnts):
        if i < 7: continue
        assert(tol_funcs.tolerable(mean_res[i,0], Aeq, tolerance))
        assert(tol_funcs.tolerable(mean_res[i,1], Beq, tolerance))

        
