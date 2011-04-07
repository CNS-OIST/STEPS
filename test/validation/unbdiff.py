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

import datetime
import steps.model as smodel
import math
import numpy
import steps.solver as solvmod
import steps.utilities.meshio as smeshio
import steps.geom as stetmesh
import steps.rng as srng
import time
import pylab

from tol_funcs import *

rng = srng.create('mt19937', 1024) 
rng.initialize(int(time.time()%4294967295)) # The max unsigned long

# Number of iterations; plotting dt; sim endtime:
NITER = 10

DT = 0.01
INT = 0.21

# Number of molecules injected in centre; diff constant; number of tets sampled:
NINJECT = 100000 	

DCST = 0.02e-9

# With good code <1% fail with a tolerance of 5% 
tolerance = 5.0/100

########################################################################

SAMPLE = 32552	 # all tets

MESHFILE = 'sphere_rad10_33Ktets_adaptive'

# create the array of tet indices to be found at random
tetidxs = numpy.zeros(SAMPLE, dtype = 'int')
for i in range(SAMPLE): tetidxs[i] = i

# further create the array of tet barycentre distance to centre
tetrads = numpy.zeros(SAMPLE)
tetvols = numpy.zeros(SAMPLE)

########################################################################

def gen_model():
   
    mdl = smodel.Model()
    X = smodel.Spec('X', mdl)
    
    cytosolv = smodel.Volsys('cytosolv', mdl)
    dif_X = smodel.Diff('diffX', cytosolv, X, dcst = DCST)

    return mdl
	
########################################################################

def gen_geom():
    mesh = smeshio.loadMesh('validation/meshes/'+MESHFILE)[0]
    ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.0])

    ntets = mesh.countTets()
    comp = stetmesh.TmComp('cyto', mesh, xrange(ntets))
    comp.addVolsys('cytosolv')
    
    # Now find the distance of the centre of the tets to the centre of the centre tet (at 0,0,0)
    cbaryc = mesh.getTetBarycenter(ctetidx)
    for i in range(SAMPLE):
        baryc = mesh.getTetBarycenter(tetidxs[i])
        r2 = math.pow((baryc[0]-cbaryc[0]),2) + math.pow((baryc[1]-cbaryc[1]),2) + math.pow((baryc[2]-cbaryc[2]),2)
        r = math.sqrt(r2)
        # Conver to microns
        tetrads[i] = r*1.0e6
        tetvols[i] = mesh.getTetVol(tetidxs[i])
    
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

# Fetch the index of the centre tet
ctetidx = g.findTetByPoint([0.0, 0.0, 0.0])
# And fetch the total number of tets to make the data structures
ntets = g.countTets()

sim = solvmod.Tetexact(m, g, rng)

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

#Create the big old data structure: iterations x time points x concentrations
res = numpy.zeros((NITER, ntpnts, SAMPLE))

for j in xrange(NITER):
    sim.reset()
    sim.setTetCount(ctetidx, 'X', NINJECT)
    for i in range(ntpnts):
        sim.run(tpnts[i])
        for k in range(SAMPLE):
            res[j, i, k] = sim.getTetCount(tetidxs[k], 'X')
    #print '%d / %d' % (j + 1, NITER)

itermeans = numpy.mean(res, axis = 0)


tpnt_compare = [10, 15, 20]
passed = True
max_err = 0.0

for t in tpnt_compare:
    bin_n = 20
    
    r_max = tetrads.max()
    r_min = 0.0
    
    r_seg = (r_max-r_min)/bin_n
    bin_mins = numpy.zeros(bin_n+1)
    r_tets_binned = numpy.zeros(bin_n)
    bin_vols = numpy.zeros(bin_n)    
    
    r = r_min
    for b in range(bin_n + 1):
        bin_mins[b] = r
        if (b!=bin_n): r_tets_binned[b] = r +r_seg/2.0
        r+=r_seg
    bin_counts = [None]*bin_n
    for i in range(bin_n): bin_counts[i] = []
    for i in range((itermeans[t].size)):
        i_r = tetrads[i]
        for b in xrange(bin_n):
            if(i_r>=bin_mins[b] and i_r<bin_mins[b+1]):
                bin_counts[b].append(itermeans[t][i])
                bin_vols[b]+=sim.getTetVol(tetidxs[i])
                break
                
    bin_concs = numpy.zeros(bin_n)
    for c in range(bin_n): 
        for d in range(bin_counts[c].__len__()):
            bin_concs[c] += bin_counts[c][d]
        bin_concs[c]/=(bin_vols[c]*1.0e18)
    
    for i in range(bin_n):
        if (r_tets_binned[i] > 2.0 and r_tets_binned[i] < 6.0):
            rad = r_tets_binned[i]*1.0e-6
            det_conc = 1e-18*((NINJECT/(math.pow((4*math.pi*DCST*tpnts[t]),1.5)))*(math.exp((-1.0*(rad*rad))/(4*DCST*tpnts[t]))))
            steps_conc = bin_concs[i]
            if not tolerable(det_conc, steps_conc, tolerance): passed = False
"""            if (abs(2*(det_conc-steps_conc)/(det_conc+steps_conc)) > max_err): max_err = abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))
    print max_err*100.0, "%"        
print "max error", max_err*100.0, "%"            
"""        
