########################################################################

# Stochastic second-order irreversible reaction: [A]0 != [B]0
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Reaction model (see validation/second_order_irev_AB.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import math
import time 
import numpy
import random

from . import tol_funcs

print("Reaction - Second order, irreversible, A0!=B0:")
from . import second_order_irev_AB_cp

########################################################################

def test_soirevAB():
    KCST = 5.0e6			# The reaction constant

    CONCA = 1.0e-6
    n = 2
    CONCB = CONCA/n
    VOL = 9.0e-18			

    NITER = 1000			# The number of iterations
    DT = 0.1			# Sampling time-step
    INT = 1.1			# Sim endtime

    # In test runs, with good code, <0.1% will fail with a tolerance of 1% 
    tolerance = 1.0/100

    ########################################################################

    mdl = Model()
    r = ReactionManager()
    with mdl:
        SA, SB, SC = Species.Create()
        volsys = VolumeSystem.Create()
        with volsys:
            SA + SB >r['R1']> SC
            r['R1'].K = KCST

    geom = Geometry()
    with geom:
        comp1 = Compartment.Create(volsys, VOL)

    rng = RNG('mt19937', 512, int(random.random()*1000))

    sim = Simulation('Wmdirect', mdl, geom, rng)

    rs = ResultSelector(sim)

    res = rs.comp1.LIST(SA, SB).Conc

    sim.toSave(res, dt=DT)

    for i in range (0, NITER):
        sim.newRun()
        sim.restore('./validation_cp/cp/second_order_irev_AB')
        sim.run(INT)

    mean_res = numpy.mean(res.data, 0)

    SC = CONCA-CONCB
    for t, (SA, SB) in zip(res.time[0], mean_res):
        lnBA = math.log(SB/SA)
        lineAB = math.log(CONCB/CONCA) -SC*KCST*t
        
        assert(tol_funcs.tolerable(lnBA, lineAB, tolerance))


