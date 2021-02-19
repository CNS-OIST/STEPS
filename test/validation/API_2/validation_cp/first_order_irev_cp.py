########################################################################

# Stochastic first-order irreversible reaction.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Irreversible 
# Reaction model (see validation/first_order_irev.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import numpy
import time

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

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res_m = numpy.zeros([NITER, ntpnts, 1])
res_std1 = numpy.zeros([ntpnts, 1])
res_std2 = numpy.zeros([ntpnts, 1])

sim.newRun()
sim.comp1.SA.Count = N
sim.checkpoint('./validation_cp/cp/first_order_irev')


