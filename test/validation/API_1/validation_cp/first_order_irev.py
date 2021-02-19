########################################################################

# Stochastic first-order irreversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Irreversible 
# Reaction model (see validation/first_order_irev.py)
  
########################################################################

from __future__ import print_function, absolute_import
import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import numpy as np
import time

from . import tol_funcs

print("Reaction - First order, irreversible:")
from . import first_order_irev_cp

########################################################################

KCST = 5			# The reaction constant
N = 50				# Can set count or conc
VOL = 1.0e-18

NITER = 100000		# The number of iterations
DT = 0.1			# Sampling time-step
INT = 1.1			# Sim endtime

# Tolerance for the comparison:
# In test runs, with good code, < 1%  will fail with a 1.5% tolerance
tolerance = 2.0/100  

########################################################################

def test_foirev():
    mdl  = smod.Model()

    A = smod.Spec('A', mdl)
    volsys = smod.Volsys('vsys',mdl)
    R1 = smod.Reac('R1', volsys, lhs = [A], rhs = [], kcst = KCST)


    geom = sgeom.Geom()
    comp1 = sgeom.Comp('comp1', geom, VOL)
    comp1.addVolsys('vsys')

    rng = srng.create('mt19937', 1000)
    rng.initialize(int(time.time()%4294967295))


    sim = ssolv.Wmdirect(mdl, geom, rng)
    sim.reset()

    tpnts = np.arange(0.0, INT, DT)
    ntpnts = tpnts.shape[0]

    res_m = np.zeros([NITER, ntpnts, 1])
    res_std1 = np.zeros([ntpnts, 1])
    res_std2 = np.zeros([ntpnts, 1])

    for i in range (0, NITER):
        sim.restore('./validation_cp/cp/first_order_irev')
        for t in range(0, ntpnts):
            sim.run(tpnts[t])
            res_m[i, t, 0] = sim.getCompCount('comp1', 'A')

    mean_res = np.mean(res_m, 0)
    std_res = np.std(res_m, 0)

    m_tol = 0
    s_tol=0

    passed = True
    for i in range(ntpnts):
        if i == 0: continue
        analy = N*np.exp(-KCST*tpnts[i])
        std = np.power((N*(np.exp(-KCST*tpnts[i]))*(1-(np.exp(-KCST*tpnts[i])))), 0.5)
        if not tol_funcs.tolerable(analy, mean_res[i], tolerance):
            passed = False
        assert(tol_funcs.tolerable(std, std_res[i], tolerance))


