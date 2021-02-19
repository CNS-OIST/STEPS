########################################################################

# Stochastic degradation-diffusion process.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the spatial stochastic 
# solver 'Tetexact' in the context of the Degradation-Diffusion model 
# (see validation/kisilevich.py)
  
########################################################################
import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import math
import time 
import numpy

########################################################################

NITER =	50			# The number of iterations
DT = 0.1			# Sampling time-step
INT = 0.3		# Sim endtime

DCSTA = 400*1e-12
DCSTB = DCSTA
RCST = 100000.0e6

#NA0 = 100000	# 1000000			# Initial number of A molecules
NA0 = 1000
NB0 = NA0		# Initial number of B molecules

SAMPLE = 1686

# <1% fail with a tolerance of 7.5%
tolerance = 7.5/100


########################################################################

mdl = Model()
r = ReactionManager()
with mdl:
    SA, SB = Species.Create()
    volsys = VolumeSystem.Create()
    with volsys:
        SA + SB >r['R1']> None
        r['R1'].K = RCST
        
        D_a =     Diffusion.Create(SA, DCSTA)
        D_b =     Diffusion.Create(SB, DCSTB)


mesh = TetMesh.Load('./validation_rd/meshes/brick_40_4_4_1686tets')

with mesh:
    acomptets = TetList(tet for tet in mesh.tets if tet.center.x < 0)
    bcomptets = mesh.tets - acomptets
    
    compa = Compartment.Create(acomptets, volsys)
    compb = Compartment.Create(bcomptets, volsys)
    
    diffb = DiffBoundary.Create(acomptets.surface & bcomptets.surface)


rng = RNG('mt19937', 16384, int(time.time()%4294967295))

sim = Simulation('Tetexact', mdl, mesh, rng)

sim.newRun()

sim.diffb.SA.DiffusionActive = True
sim.diffb.SB.DiffusionActive = True

sim.compa.SA.Count = NA0
sim.compb.SB.Count = NB0

sim.checkpoint('./validation_cp/cp/kisilevich')

