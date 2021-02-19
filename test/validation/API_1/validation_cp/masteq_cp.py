########################################################################

# Stochastic production and degradation well-mixed reactions.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Production-Degradation model 
# (see validation/masteq.py)
  
########################################################################

import math
import numpy
import time 

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

from . import tol_funcs

########################################################################

KCST_f = 100/(6.022e23*1.0e-15) 			# The reaction constant, production
KCST_b = 10			# The reaction constant, degradation
VOL = 1.0e-18

DT = 0.1			# Sampling time-step
INT = 200000.1 		# Sim endtime

# Tolerance for the comparison:
# In tests with good code <1% fail with tolerance of 1.5%
tolerance = 1.5/100  

########################################################################

mdl  = smod.Model()

A = smod.Spec('A', mdl)

volsys = smod.Volsys('vsys',mdl)

# Production
R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A], kcst = KCST_f)
R2 = smod.Reac('R2', volsys, lhs = [A], rhs = [], kcst = KCST_b)

geom = sgeom.Geom()

comp1 = sgeom.Comp('comp1', geom, VOL)
comp1.addVolsys('vsys')

rng = srng.create('mt19937', 1000)
rng.initialize(int(time.time()%4294967295))

sim = ssolv.Wmdirect(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res = numpy.zeros([ntpnts])

sim.reset()
sim.setCompCount('comp1', 'A', 0)

sim.checkpoint('./validation_cp/cp/masteq')

