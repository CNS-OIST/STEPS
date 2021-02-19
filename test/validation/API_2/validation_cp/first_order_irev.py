########################################################################

# Stochastic first-order irreversible reaction.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Irreversible 
# Reaction model (see validation/first_order_irev.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

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
    mdl = Model()
    r = ReactionManager()
    with mdl:
        SA = Species.Create()
        volsys = VolumeSystem.Create()
        with volsys:
            SA >r['R1']> None
            r['R1'].K = KCST

    geom = Geometry()
    with geom:
        comp1 = Compartment.Create(volsys, VOL)

    rng = RNG('mt19937', 1000, int(time.time()%4294967295))

    sim = Simulation('Wmdirect', mdl, geom, rng)

    rs = ResultSelector(sim)

    res = rs.comp1.SA.Count

    sim.toSave(res, dt=DT)

    for i in range (0, NITER):
        sim.newRun()
        sim.restore('./validation_cp/cp/first_order_irev')
        sim.run(INT)

    mean_res = np.mean(res.data, 0)
    std_res = np.std(res.data, 0)

    m_tol = 0
    s_tol=0

    passed = True
    for i in range(len(res.time[0])):
        if i == 0:
            continue
        analy = N*np.exp(-KCST*res.time[0,i])
        std = np.power((N*(np.exp(-KCST*res.time[0,i]))*(1-(np.exp(-KCST*res.time[0,i])))), 0.5)
        if not tol_funcs.tolerable(analy, mean_res[i], tolerance):
            passed = False
        assert(tol_funcs.tolerable(std, std_res[i], tolerance))


