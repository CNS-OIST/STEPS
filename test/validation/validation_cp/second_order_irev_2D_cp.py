########################################################################

# Stochastic second-order irreversible reaction on a surface.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Surface Reaction model (see validation/second_order_irev_2D.py)
  
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

VOL = 1.0e-18			

COUNTA = 100.0
n=2.0
COUNTB = COUNTA/n 


KCST = 10.0e10			# The reaction constant

AREA = 10.0e-12

CCST = KCST/(6.02214179e23*AREA)


NITER = 1000			# The number of iterations
DT = 0.05			# Sampling time-step
INT = 1.05			# Sim endtime

# In tests fewer than 0.1% fail with tolerance of 2%
tolerance = 2.0/100

########################################################################

mdl  = smod.Model()

A = smod.Spec('A', mdl)
B = smod.Spec('B', mdl)
C = smod.Spec('C', mdl)

surfsys = smod.Surfsys('ssys',mdl)

SR1 = smod.SReac('SR1', surfsys, slhs = [A, B], srhs = [C], kcst = KCST)

geom = sgeom.Geom()

comp1 = sgeom.Comp('comp1', geom, VOL)
patch1 = sgeom.Patch('patch1', geom, comp1, area = AREA)
patch1.addSurfsys('ssys')

import random
rng = srng.create('mt19937', 1000)
rng.initialize(int(random.random()*4294967295))


sim = ssolv.Wmdirect(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res_m = numpy.zeros([NITER, ntpnts, 3])

sim.reset()
sim.setPatchCount('patch1', 'A', COUNTA)
sim.setPatchCount('patch1', 'B', COUNTB)
sim.checkpoint('./validation_cp/cp/second_order_irev_2D')

