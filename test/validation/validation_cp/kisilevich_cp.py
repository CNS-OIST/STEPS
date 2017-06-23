########################################################################

# Stochastic degradation-diffusion process.
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the spatial stochastic 
# solver 'Tetexact' in the context of the Degradation-Diffusion model 
# (see validation/kisilevich.py)
  
########################################################################

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import math
import time 
import numpy
import steps.utilities.meshio as meshio

from tol_funcs import *

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


# create the array of tet indices to be found at random
tetidxs = numpy.zeros(SAMPLE, dtype = 'int')
# further create the array of tet barycentre distance to centre
tetrads = numpy.zeros(SAMPLE)

########################################################################
rng = srng.create('mt19937', 512) 
rng.initialize(int(time.time()%4294967295)) # The max unsigned long


mdl  = smod.Model()

A = smod.Spec('A', mdl)
B = smod.Spec('B', mdl)

volsys = smod.Volsys('vsys',mdl)


R1 = smod.Reac('R1', volsys, lhs = [A,B], rhs = [])

R1.setKcst(RCST)

D_a = smod.Diff('D_a', volsys, A)
D_a.setDcst(DCSTA)
D_b = smod.Diff('D_b', volsys, B)
D_b.setDcst(DCSTB)


mesh = meshio.loadMesh('./validation_rd/meshes/brick_40_4_4_1686tets')[0]

VOLA = mesh.getMeshVolume()/2.0
VOLB = VOLA

ntets = mesh.countTets()

acomptets = []
bcomptets = []
max = mesh.getBoundMax()
min = mesh.getBoundMax()
midz = 0.0
compatris=set()
compbtris=set()
for t in range(ntets):
    barycz = mesh.getTetBarycenter(t)[0]
    tris = mesh.getTetTriNeighb(t)
    if barycz < midz: 
        acomptets.append(t)
        compatris.add(tris[0])
        compatris.add(tris[1])
        compatris.add(tris[2])
        compatris.add(tris[3])
    else: 
        bcomptets.append(t)
        compbtris.add(tris[0])
        compbtris.add(tris[1])
        compbtris.add(tris[2])
        compbtris.add(tris[3])

dbset = compatris.intersection(compbtris)
dbtris = list(dbset)

compa = sgeom.TmComp('compa', mesh, acomptets)
compb = sgeom.TmComp('compb', mesh, bcomptets)
compa.addVolsys('vsys')
compb.addVolsys('vsys')

diffb = sgeom.DiffBoundary('diffb', mesh, dbtris)


# Now fill the array holding the tet indices to sample at random
assert(SAMPLE <= ntets)

numfilled = 0
while (numfilled < SAMPLE):
    tetidxs[numfilled] = numfilled
    numfilled +=1

# Now find the distance of the centre of the tets to the Z lower face
for i in range(SAMPLE):
	baryc = mesh.getTetBarycenter(int(tetidxs[i]))
	r = baryc[0]
	tetrads[i] = r*1.0e6

Atets = acomptets
Btets = bcomptets

rng = srng.create('mt19937', 16384)
rng.initialize(int(time.time()%4294967295))


sim = ssolv.Tetexact(mdl, mesh, rng)

sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

resA = numpy.zeros((NITER, ntpnts, SAMPLE))
resB = numpy.zeros((NITER, ntpnts, SAMPLE))

sim.reset()

sim.setDiffBoundaryDiffusionActive('diffb', 'A', True)
sim.setDiffBoundaryDiffusionActive('diffb', 'B', True)

sim.setCompCount('compa', 'A', NA0)
sim.setCompCount('compb', 'B', NB0)

sim.checkpoint('./validation_cp/cp/kisilevich')

