########################################################################

# Stochastic first-order reversible reaction.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Reversible 
# Reaction model (see validation/first_order_rev.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import time 
import numpy

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

sim.newRun()
sim.comp1.SA.Count = COUNT
sim.comp1.SB.Count = 0.0
sim.checkpoint('./validation_cp/cp/first_order_rev')



