########################################################################

# Stochastic first-order reversible reaction.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Reversible 
# Reaction model (see validation/first_order_rev.py)
  
########################################################################

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import time 
import numpy

from . import tol_funcs

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

sim.reset()
sim.setCompCount('comp1', 'A', COUNT)
sim.setCompCount('comp1', 'B', 0.0)
sim.checkpoint('./validation_cp/cp/first_order_rev')



