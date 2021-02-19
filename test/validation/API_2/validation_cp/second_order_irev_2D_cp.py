########################################################################

# Stochastic second-order irreversible reaction on a surface.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Surface Reaction model (see validation/second_order_irev_2D.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import math
import numpy
import random
import time 

from scipy.constants import Avogadro

########################################################################

VOL = 1.0e-18			

COUNTA = 100.0
n=2.0
COUNTB = COUNTA/n 


KCST = 10.0e10			# The reaction constant

AREA = 10.0e-12

CCST = KCST/(Avogadro*AREA)


NITER = 1000			# The number of iterations
DT = 0.05			# Sampling time-step
INT = 1.05			# Sim endtime

# In tests fewer than 0.1% fail with tolerance of 2%
tolerance = 2.0/100

########################################################################

mdl = Model()
r = ReactionManager()
with mdl:
    SA, SB, SC = Species.Create()
    surfsys = SurfaceSystem.Create()
    with surfsys:
        SA.s + SB.s >r['SR1']> SC.s
        r['SR1'].K = KCST

geom = Geometry()
with geom:
    comp1 = Compartment.Create(None, VOL)
    patch1 = Patch.Create(comp1, None, surfsys, AREA)

rng = RNG('mt19937', 1000, int(random.random()*4294967295))

sim = Simulation('Wmdirect', mdl, geom, rng)

sim.newRun()
sim.patch1.SA.Count = COUNTA
sim.patch1.SB.Count = COUNTB
sim.checkpoint('./validation_cp/cp/second_order_irev_2D')

