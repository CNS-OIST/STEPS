# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2011 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import math
import time 
import numpy
from pylab import *
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

D_a = smod.Diff('D_a', volsys, A, dcst = DCSTA)
D_b = smod.Diff('D_b', volsys, B, dcst = DCSTB)


mesh = meshio.loadMesh('validation/meshes/brick_40_4_4_1686tets')[0]

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
	baryc = mesh.getTetBarycenter(tetidxs[i])
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


for i in xrange (0, NITER):
    sim.reset()
    
    sim.setDiffBoundaryDiffusionActive('diffb', 'A', True)
    sim.setDiffBoundaryDiffusionActive('diffb', 'B', True)

    sim.setCompCount('compa', 'A', NA0)
    sim.setCompCount('compb', 'B', NB0)

    for t in xrange(0, ntpnts):
        sim.run(tpnts[t])
        for k in range(SAMPLE):
            resA[i,t,k] = sim.getTetCount(tetidxs[k], 'A')
            resB[i,t,k] = sim.getTetCount(tetidxs[k], 'B')
	
	
itermeansA = numpy.mean(resA, axis=0)
itermeansB = numpy.mean(resB, axis=0)

import pylab


def getdetc(t, x):
    N = 1000		# The number to represent infinity in the exponential calculation
    L = 20e-6

    concA  = 0.0
    for n in xrange(N):
        concA+= ((1.0/(2*n +1))*math.exp((-(DCSTA/(20.0e-6))*math.pow((2*n +1), 2)*math.pow(math.pi, 2)*t)/(4*L))*math.sin(((2*n +1)*math.pi*x)/(2*L)))
    concA*=((4*NA0/math.pi)/(VOLA*6.022e26))*1.0e6	
    
    return concA


tpnt_compare = [1, 2]
passed = True
max_err = 0.0

for tidx in tpnt_compare:
    NBINS=10
    radmax = 0.0
    radmin = 10.0
    for r in tetrads:
        if (r > radmax): radmax = r
        if (r < radmin) : radmin = r
    
    rsec = (radmax-radmin)/NBINS
    binmins = numpy.zeros(NBINS+1)
    tetradsbinned = numpy.zeros(NBINS)
    r = radmin
    bin_vols = numpy.zeros(NBINS)
    
    for b in xrange(NBINS+1):
        binmins[b] = r
        if (b!=NBINS): tetradsbinned[b] = r +rsec/2.0
        r+=rsec
    
    bin_countsA = [None]*NBINS
    bin_countsB = [None]*NBINS
    for i in xrange(NBINS):
        bin_countsA[i] = []
        bin_countsB[i] = []
    filled = 0
    
    for i in xrange(itermeansA[tidx].size):
        irad = tetrads[i]
        
        for b in xrange(NBINS):
            if(irad>=binmins[b] and irad<binmins[b+1]):
                bin_countsA[b].append(itermeansA[tidx][i])
                bin_vols[b]+=sim.getTetVol(tetidxs[i])
                filled+=1.0
                break
    filled = 0
    for i in xrange(itermeansB[tidx].size):
        irad = tetrads[i]
        
        for b in xrange(NBINS):
            if(irad>=binmins[b] and irad<binmins[b+1]):
                bin_countsB[b].append(itermeansB[tidx][i])
                filled+=1.0
                break
    
    bin_concsA = numpy.zeros(NBINS)
    bin_concsB = numpy.zeros(NBINS)
    
    for c in range(NBINS): 
        for d in range(bin_countsA[c].__len__()):
            bin_concsA[c] += bin_countsA[c][d]
        for d in range(bin_countsB[c].__len__()):
            bin_concsB[c] += bin_countsB[c][d]        
        
        bin_concsA[c]/=(bin_vols[c])
        bin_concsA[c]*=(1.0e-3/6.022e23)*1.0e6    
        bin_concsB[c]/=(bin_vols[c])
        bin_concsB[c]*=(1.0e-3/6.022e23)*1.0e6 
    
    for i in range(NBINS):
        rad = abs(tetradsbinned[i])*1.0e-6
        
        if (tetradsbinned[i] < -5):
            # compare A
            det_conc = getdetc(tpnts[tidx], rad)
            steps_conc = bin_concsA[i]
            if not tolerable(det_conc, steps_conc, tolerance): passed = False
            if (abs(2*(det_conc-steps_conc)/(det_conc+steps_conc)) > max_err): max_err = abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))
            #print "Error:",abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))*100.0, "%"    
            #print ""    
        elif (tetradsbinned[i] > 5):
            # compare B
            det_conc = getdetc(tpnts[tidx], rad)
            steps_conc = bin_concsB[i]
            if not tolerable(det_conc, steps_conc, tolerance): passed = False
            if (abs(2*(det_conc-steps_conc)/(det_conc+steps_conc)) > max_err): max_err = abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))
            #print "Error:",abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))*100.0, "%"       
            #print "" 
#print "max error", max_err*100.0, "%"               
