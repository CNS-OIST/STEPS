########################################################################

# Stochastic first-order irreversible reaction.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the First-Order Irreversible 
# Reaction model (see validation/first_order_irev.py)
  
########################################################################

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import numpy
import time

from tol_funcs import *

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

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res_m = numpy.zeros([NITER, ntpnts, 1])
res_std1 = numpy.zeros([ntpnts, 1])
res_std2 = numpy.zeros([ntpnts, 1])

sim.reset()
sim.setCompCount('comp1', 'A', N)
sim.checkpoint('./validation_cp/cp/first_order_irev')


