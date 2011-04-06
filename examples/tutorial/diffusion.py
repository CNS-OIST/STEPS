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

# Example: Unbounded diffusion
# http://steps.sourceforge.net/manual/diffusion.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import numpy
import pylab
import random
import time

import steps.model as smodel
import steps.solver as solvmod
import steps.geom as stetmesh
import steps.rng as srng

import steps.utilities.meshio as smeshio

########################################################################

# The number of iterations to run 
NITER = 10
# The data collection time increment (s)
DT = 0.001
# The simulation endtime (s)
INT = 0.101

# The number of molecules to be injected into the centre
NINJECT = 10000	

# The number of tetrahedral elements to sample data from. 
SAMPLE = 2000	

# The diffusion constant for our diffusing species (m^2/s)
DCST= 20.0e-12

# Array to hold tetrahedron indices (integers)
tetidxs = numpy.zeros(SAMPLE, dtype = 'int')

# Array to hold tetrahedron radial distances from mesh centre (floats) 
tetrads = numpy.zeros(SAMPLE)

beg_time = time.time()

########################################################################

def printtime(end_time):

	totsecs = int(end_time-beg_time)
	sec = totsecs%60
	totmin = totsecs/60
	min = totmin%60
	hours = totmin/60
	
	print 'Simulation time: %d h, %d min, %d sec' %(hours, min, sec)
	
########################################################################

def gen_model():         
	mdl = smodel.Model()        
	A = smodel.Spec('A', mdl)                        
	vsys = smodel.Volsys('cytosolv', mdl)
	diff_A = smodel.Diff('diff_A', vsys, A, dcst = DCST)
	return mdl
	
########################################################################

def gen_geom():
    
    print "Loading mesh..."
    mesh = smeshio.loadMesh('meshes/sphere_rad10_11Ktets')[0]
    print "Mesh Loaded"
    
    # Find the total number of tetrahedrons in the mesh	
    ntets = mesh.countTets()
    # Create a compartment containing all tetrahedron
    comp = stetmesh.TmComp('cyto', mesh, range(ntets))
    comp.addVolsys('cytosolv')
    
    print "Finding tetrahedron samples..."
    # Fetch the central tetrahedron index and store:
    ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.0])
    tetidxs[0] = ctetidx

    # Find the central tetrahedron's four neighbours:
    neighbidcs = mesh.getTetTetNeighb(ctetidx)
    tetidxs[1], tetidxs[2], tetidxs[3], tetidxs[4] = neighbidcs
    
    # Keep track how many tet indices we have stored so far
    stored = 5
    
    # Find the maximum and minimum coordinates of the mesh
    max = mesh.getBoundMax()
    min = mesh.getBoundMin()
    
    # Run a loop until we have stored all tet indices we require
    while (stored < SAMPLE):

        # Fetch 3 random numbers between 0 and 1
        rnx = random.random()
        rny = random.random()
        rnz = random.random()  
              
        # Find the coordinates in the mesh that these numbers relate to
        xcrd = min[0] + (max[0]-min[0])*rnx
        ycrd = min[1] + (max[1]-min[1])*rny
        zcrd = min[2] + (max[2]-min[2])*rnz
        
        # Find the tetrahedron that encompasses this point.
        tidx = mesh.findTetByPoint([xcrd, ycrd, zcrd])
        
        # -1 was returned if point is outside the mesh:
        if (tidx == -1): continue
        if (tidx not in tetidxs):
            tetidxs[stored] = tidx
            stored += 1
    
    # Find the barycenter of the central tetrahedron
    cbaryc = mesh.getTetBarycenter(ctetidx)
    
    for i in range(SAMPLE):
        # Fetch the barycenter of the tetrahedron:
        baryc = mesh.getTetBarycenter(tetidxs[i])

        # Find the radial distance of this tetrahedron to mesh center:
        r = math.sqrt(math.pow((baryc[0]-cbaryc[0]),2) + \
                math.pow((baryc[1]-cbaryc[1]),2) + \
                    math.pow((baryc[2]-cbaryc[2]),2))

        # Store the radial distance (in microns):
        tetrads[i] = r*1.0e6
        
    print "Tetrahedron samples found"
    
    return mesh
    
########################################################################

model = gen_model()
tmgeom = gen_geom()

rng = srng.create('mt19937', 512) 
rng.initialize(2903) # The max unsigned long

sim = solvmod.Tetexact(model, tmgeom, rng)

tpnts = numpy.arange(0.0, INT, DT)
# Find how many "time points" we have
ntpnts = tpnts.shape[0]

# Create the data structure: iterations x time points x tet samples
res = numpy.zeros((NITER, ntpnts, SAMPLE))

# Fetch the index of the tetrahedron at the centre of the mesh
ctetidx = tmgeom.findTetByPoint([0.0, 0.0, 0.0])

# Run NITER number of iterations:
for i in range(NITER):
	sim.reset()
	print "Running iteration", i
    # Inject all molecules into the central tet:
	sim.setTetCount(ctetidx, 'A', NINJECT)
	for j in range(ntpnts):
		sim.run(tpnts[j])
        # Loop over the tetrahedrons we are saving data for
		for k in range(SAMPLE):
            # Save the concentration in the tetrahedron, in uM
			res[i, j, k] = sim.getTetConc(tetidxs[k], 'A')*1.0e6
	printtime(time.time())

res_mean = numpy.mean(res, axis = 0)

########################################################################

def plotres(tidx):
	if (tidx >= INT/DT):
		print "Time index is out of range."
		return
	
	pylab.scatter(tetrads, res_mean[tidx], s=2)
	pylab.xlabel('Radial distance of tetrahedron ($\mu$m)')			
	pylab.ylabel('Concentration in tetrahedron ($\mu$M)')
	t = tpnts[tidx]
	pylab.title('Unbounded diffusion. Time: ' + str(t) + 's')
	plotanlyt(t)
	pylab.xlim(0.0, 10.0)
	pylab.ylim(0.0)
	pylab.show()

########################################################################

def plotanlyt(t): 	
	segs = 100 	
	anlytconc = numpy.zeros((segs)) 	
	radialds = numpy.zeros((segs)) 	
	maxrad = 0.0 	
	for i in tetrads: 		
		if (i > maxrad): maxrad = i 	
	maxrad *= 1e-6 	
	intervals = maxrad/segs 	
	rad = 0.0 	
	for i in range((segs)): 		
		# Find the conc from analytical solution, and convert to mol/L 		
		anlytconc[i]=1.0e3*(1/6.022e23)* \
            ((NINJECT/(math.pow((4*math.pi*DCST*t),1.5)))* \
                (math.exp((-1.0*(rad*rad))/(4*DCST*t)))) 		
		radialds[i] = rad*1e6 		
		rad += intervals 	
	pylab.plot(radialds, anlytconc, color = 'red')
	
########################################################################

plotres(100)
