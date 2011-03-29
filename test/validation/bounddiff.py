import datetime
import steps.model as smodel
import math
import steps.solver as solvmod
import steps.utilities.meshio as meshio
import steps.geom as stetmesh
import steps.rng as srng
import time
import numpy
import pylab

from tol_funcs import *

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
    dif_X = smodel.Diff('diffX', cytosolv, X, dcst = DCST)
    
    return mdl
	
########################################################################

def gen_geom():
    mesh = meshio.loadMesh('validation/meshes/' +MESHFILE)[0]
    
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
        baryc = mesh.getTetBarycenter(tetidxs[i])
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

for i in xrange(ntets):
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
for i in xrange(boundtets.__len__()):
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

for j in xrange(NITER):
    sim.reset()
    tetcount = int((1.0*NINJECT)/nztets)
    totset = 0
    for k in minztets:
        sim.setTetCount(k, 'X', tetcount)
        totset+=tetcount
    for i in range(ntpnts):
        sim.run(tpnts[i])
        for k in range(SAMPLE):
            res[j, i, k] = sim.getTetCount(tetidxs[k], 'X')
    #print '%d / %d' % (j + 1, NITER)

itermeans = numpy.mean(res, axis = 0)
	

########################################################################

D = DCST
pi = math.pi
nmax = 1000
N = NINJECT
N = int((1.0*NINJECT)/nztets)*nztets
def getprob(x,t):
	if(x>a): 
		print 'x out of bounds'
		return
	p=0.0
	for n in xrange(nmax):
		if (n==0): A = math.sqrt(1.0/a)
		else : A = math.sqrt(2.0/a)
		p+= math.exp(-D*math.pow((n*pi/a), 2)*t)*A*math.cos(n*pi*x/a)*A*a
	
	return p*N/a

tpnt_compare = [6, 8, 10]
passed = True
max_err = 0.0

for t in tpnt_compare:
    NBINS = 5
    
    radmax = 0.0
    radmin = 11.0
    for r in tetrads:
        if (r > radmax): radmax = r
        if (r < radmin) : radmin = r

    rsec = (radmax-radmin)/NBINS
    binmins = numpy.zeros(NBINS+1)
    tetradsbinned = numpy.zeros(NBINS)
    r = radmin
    bin_vols = numpy.zeros(NBINS)

    for b in range(NBINS+1):
        binmins[b] = r
        if (b!=NBINS): tetradsbinned[b] = r +rsec/2.0
        r+=rsec
    
    bin_counts = [None]*NBINS
    for i in xrange(NBINS):
        bin_counts[i] = []
    filled = 0

    for i in range(itermeans[t].size):
        irad = tetrads[i]
        
        for b in xrange(NBINS):
            if(irad>=binmins[b] and irad<binmins[b+1]):
                bin_counts[b].append(itermeans[t][i])
                bin_vols[b]+=sim.getTetVol(tetidxs[i])
                filled+=1.0
                break
    bin_concs = numpy.zeros(NBINS)
    for c in range(NBINS): 
        for d in range(bin_counts[c].__len__()):
            bin_concs[c] += bin_counts[c][d]
        bin_concs[c]/=(bin_vols[c])
        bin_concs[c]*=(1.0e-3/6.022e23)*1.0e6
    
    for i in range(NBINS):
        if (tetradsbinned[i] > 2 and tetradsbinned[i] < 8):
            rad = tetradsbinned[i]*1.0e-6
            det_conc = (getprob(rad, tpnts[t])/area)*(1.0/6.022e20)
            steps_conc = bin_concs[i]
            if not tolerable(det_conc, steps_conc, tolerance): passed = False
            if (abs(2*(det_conc-steps_conc)/(det_conc+steps_conc)) > max_err): max_err = abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))
            #print "Error:",abs(2*(det_conc-steps_conc)/(det_conc+steps_conc))*100.0, "%"        
#print "max error", max_err*100.0, "%"               

