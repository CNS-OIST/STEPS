########################################################################

# 1D diffusion in a finite tube with fully-reflective boundary.
# Circular "point source" from one face. 
# CHECKPOINT

# AIMS: to verify checkpointing and restoring of the spatial stochastic 
# solver 'Tetexact' in the context of the Bounded Diffusion model 
# (see validation/bounddiff.py)
  
########################################################################

import datetime
import steps.model as smodel
import math
import steps.solver as solvmod
import steps.utilities.meshio as meshio
import steps.geom as stetmesh
import steps.rng as srng
import time
import numpy

from . import tol_funcs

rng = srng.create('mt19937', 512) 
rng.initialize(int(time.time()%4294967295)) # The max unsigned long

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

# create the array of tet indices to be found at random
tetidxs = numpy.zeros(SAMPLE, dtype = 'int')
# further create the array of tet barycentre distance to centre
tetrads = numpy.zeros(SAMPLE)


########################################################################

def gen_model():
    
    mdl = smodel.Model()
    
    X = smodel.Spec('X', mdl)
    cytosolv = smodel.Volsys('cytosolv', mdl)
    dif_X = smodel.Diff('diffX', cytosolv, X)
    dif_X.setDcst(DCST)
    
    return mdl

########################################################################

def gen_geom():
    mesh = meshio.loadMesh('./validation_rd/meshes/' +MESHFILE)[0]
    
    a = mesh.getBoundMax()[2]-mesh.getBoundMin()[2]
    area = mesh.getMeshVolume()/a
    
    ntets = mesh.countTets()
    comp = stetmesh.TmComp('cyto', mesh, range(ntets))
    comp.addVolsys('cytosolv')
    
    assert(SAMPLE == ntets)
    
    numfilled = 0
    while (numfilled < SAMPLE):
        tetidxs[numfilled] = numfilled
        numfilled +=1
    
    # Now find the distance of the centre of the tets to the Z lower face
    for i in range(SAMPLE):
        baryc = mesh.getTetBarycenter(int(tetidxs[i]))
        min = mesh.getBoundMin()
        r = baryc[2] - min[2]
        # Convert to microns
        tetrads[i] = r*1.0e6
    
    return mesh, area, a

########################################################################

m = gen_model()
g, area, a = gen_geom()

# And fetch the total number of tets to make the data structures
ntets = g.countTets()

sim = solvmod.Tetexact(m, g, rng)

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]


# Create the big old data structure: iterations x time points x concentrations
res = numpy.zeros((NITER, ntpnts, SAMPLE))

# Find the tets connected to the bottom face
# First find all the tets with ONE face on a boundary
boundtets = []
#store the 0to3 index of the surface triangle for each of these boundary tets
bt_srftriidx = []

for i in range(ntets):
	tettemp = g.getTetTetNeighb(i)
	if (tettemp[0] ==-1 or tettemp[1] == -1 or tettemp[2] == -1 or tettemp[3] == -1): 
		boundtets.append(i)
		templist = []
		if (tettemp[0] == -1): 
			templist.append(0)
		if (tettemp[1] == -1): 
			templist.append(1)
		if (tettemp[2] == -1): 
			templist.append(2)
		if (tettemp[3] == -1): 
			templist.append(3)
		bt_srftriidx.append(templist)

assert (boundtets.__len__() == bt_srftriidx.__len__())

minztets = []
boundminz = g.getBoundMin()[2] + 0.01e-06
num2s=0
for i in range(boundtets.__len__()):
	# get the boundary triangle
	if (bt_srftriidx[i].__len__() == 2): num2s+=1
	for btriidx in bt_srftriidx[i]:
		zminboundtri = True
		tribidx = g.getTetTriNeighb(boundtets[i])[btriidx]
		tritemp = g.getTri(tribidx)
		trizs = [0.0, 0.0, 0.0]
		trizs[0] = g.getVertex(tritemp[0])[2]
		trizs[1] = g.getVertex(tritemp[1])[2]
		trizs[2] = g.getVertex(tritemp[2])[2]
		for j in range(3):
			if (trizs[j]>boundminz): zminboundtri = False
		if (zminboundtri): minztets.append(boundtets[i])

nztets = minztets.__len__()
volztets = 0.0
for z in minztets:
	volztets += g.getTetVol(z)
conc = NITER*6.022e23*1.0e-3/volztets

sim.reset()
tetcount = int((1.0*NINJECT)/nztets)
for k in minztets:
    sim.setTetCount(k, 'X', tetcount)
    
sim.checkpoint('./validation_cp/cp/boundiff')


