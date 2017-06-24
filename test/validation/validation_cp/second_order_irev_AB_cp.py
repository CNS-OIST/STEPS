########################################################################

# Stochastic second-order irreversible reaction: [A]0 != [B]0
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Reaction model (see validation/second_order_irev_AB.py)
  
########################################################################

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import math
import time 
import numpy

from tol_funcs import *

########################################################################

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

mdl  = smod.Model()

A = smod.Spec('A', mdl)
B = smod.Spec('B', mdl)
C = smod.Spec('C', mdl)

volsys = smod.Volsys('vsys',mdl)

R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [C], kcst = KCST)

geom = sgeom.Geom()

comp1 = sgeom.Comp('comp1', geom, VOL)
comp1.addVolsys('vsys')

rng = srng.create('mt19937', 512)
import random
#rng.initialize(int(time.time()%4294967295))
import random
rng.initialize(int(random.random()*1000))

sim = ssolv.Wmdirect(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res_m = numpy.zeros([NITER, ntpnts, 3])

sim.reset()
sim.setCompConc('comp1', 'A', CONCA)
sim.setCompConc('comp1', 'B', CONCB)
sim.checkpoint('./validation_cp/cp/second_order_irev_AB')

