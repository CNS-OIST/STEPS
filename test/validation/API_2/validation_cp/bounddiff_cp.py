########################################################################

# 1D diffusion in a finite tube with fully-reflective boundary.
# Circular "point source" from one face. 
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the spatial stochastic 
# solver 'Tetexact' in the context of the Bounded Diffusion model 
# (see validation/bounddiff.py)
  
########################################################################

import steps.interface

from steps.model import *
from steps.sim import *
from steps.geom import *
from steps.rng import *

import time
import numpy

rng = RNG('mt19937', 512, int(time.time()%4294967295))

NITER = 10
DT = 0.01
INT = 0.11

# The number of initial molecules:
NINJECT = 10000	

DCST = 0.2e-9

# In tests, with good code, <1% fail with a tolerance of 5%
tolerance = 5.0/100

# The number of tets to sample at random:
SAMPLE = 1060	

MESHFILE = 'cyl_diam2__len10_1060tets'

########################################################################

mdl = Model()
with mdl:
    X = Species.Create()
    cytosolv = VolumeSystem.Create()
    with cytosolv:
        dif_X = Diffusion.Create(X, DCST)

########################################################################

mesh = TetMesh.Load('./validation_rd/meshes/' + MESHFILE)

a = mesh.bbox.max.z - mesh.bbox.min.z
area = mesh.Vol/a

ntets = len(mesh.tets)
with mesh:
    comp = Compartment.Create(mesh.tets, mdl.cytosolv)

assert(SAMPLE == ntets)
    
########################################################################

sim = Simulation('Tetexact', mdl, mesh, rng)

boundminz = mesh.bbox.min.z + 0.01e-06
minztets = TetList(mesh=mesh)
for tri in mesh.surface:
    if all(v.z <= boundminz for v in tri.verts):
        minztets |= tri.tetNeighbs

sim.newRun()

sim.TETS(minztets).X.Count = int(NINJECT/len(minztets))
    
sim.checkpoint('./validation_cp/cp/boundiff')


