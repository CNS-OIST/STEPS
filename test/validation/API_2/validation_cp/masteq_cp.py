########################################################################

# Stochastic production and degradation well-mixed reactions.
# CHECKPOINT

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

from scipy.constants import Avogadro

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

sim.newRun()
sim.comp1.SA.Count = 0

sim.checkpoint('./validation_cp/cp/masteq')

