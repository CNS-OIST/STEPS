########################################################################

# 1D diffusion in a finite tube with fully-reflective boundary.
# Circular "point source" from one face. 
# RESTORE

# AIMS: to verify checkpointing and restoring of the spatial stochastic 
# solver 'Tetvesicle' in the context of the Bounded Diffusion model
# (see validation/bounddiff.py)
  
########################################################################

import unittest
import os
import time
import numpy
import math
from mpi4py import MPI

from steps.geom import UNKNOWN_TET
import steps.model as smodel
import steps.mpi
import steps.mpi.solver as ssolv
import steps.utilities.meshio as meshio
import steps.geom as stetmesh
import steps.rng as srng

from . import tol_funcs

comm_world = MPI.COMM_WORLD

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
# further create the array of tet barycenter distance to center
tetrads = numpy.zeros(SAMPLE)


########################################################################

class TestBDiffTetvesicleCP(unittest.TestCase):
    def setUp(self):
        myrank = steps.mpi.rank

        mdl = smodel.Model()
        
        X = smodel.Spec('X', mdl)
        cytosolv = smodel.Volsys('cytosolv', mdl)
        dif_X = smodel.Diff('diffX', cytosolv, X)
        dif_X.setDcst(DCST)

        ########################################################################

        mesh = meshio.loadMesh('./validation_rd/meshes/' +MESHFILE)[0]
        
        a = mesh.getBoundMax()[2]-mesh.getBoundMin()[2]
        area = mesh.getMeshVolume()/a
        
        ntets = mesh.countTets()
        comp = stetmesh.TmComp('cyto', mesh, range(ntets))
        comp.addVolsys('cytosolv')
        
        self.assertTrue(SAMPLE == ntets)
        
        numfilled = 0
        while (numfilled < SAMPLE):
            tetidxs[numfilled] = numfilled
            numfilled +=1
        
        # Now find the distance of the center of the tets to the Z lower face
        for i in range(SAMPLE):
            baryc = mesh.getTetBarycenter(int(tetidxs[i]))
            min = mesh.getBoundMin()
            r = baryc[2] - min[2]
            # Convert to microns
            tetrads[i] = r*1.0e6
        
        ########################################################################

        # And fetch the total number of tets to make the data structures
        ntets = mesh.countTets()

        rng = srng.create('mt19937', 512) 
        rng.initialize(int(time.time()%4294967295)) # The max unsigned long
        
        sim = ssolv.TetVesicle(mdl, mesh, rng, ssolv.EF_NONE)
        sim.setOutputSync(True, 0);
        sim.reset()
        
        sim.setVesicleDT(DT)
        
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
            tettemp = mesh.getTetTetNeighb(i)
            templist = [t for t in range(4) if tettemp[t] == UNKNOWN_TET]
            if templist:
                boundtets.append(i)
                bt_srftriidx.append(templist)

        self.assertTrue(len(boundtets) == len(bt_srftriidx))

        minztets = []
        boundminz = mesh.getBoundMin()[2] + 0.01e-06
        num2s=0
        for i in range(boundtets.__len__()):
                # get the boundary triangle
                if (bt_srftriidx[i].__len__() == 2): num2s+=1
                for btriidx in bt_srftriidx[i]:
                        zminboundtri = True
                        tribidx = mesh.getTetTriNeighb(boundtets[i])[btriidx]
                        tritemp = mesh.getTri(tribidx)
                        trizs = [0.0, 0.0, 0.0]
                        trizs[0] = mesh.getVertex(tritemp[0])[2]
                        trizs[1] = mesh.getVertex(tritemp[1])[2]
                        trizs[2] = mesh.getVertex(tritemp[2])[2]
                        for j in range(3):
                                if (trizs[j]>boundminz): zminboundtri = False
                        if (zminboundtri): minztets.append(boundtets[i])

        nztets = minztets.__len__()
        volztets = 0.0
        for z in minztets:
                volztets += mesh.getTetVol(z)
        conc = NITER*6.022e23*1.0e-3/volztets

        tetcount = int((1.0*NINJECT)/nztets)
        for k in minztets:
            sim.setTetSpecCount(k, 'X', tetcount)

        if (myrank==0):
            new_dir = './validation_cp_mpi/cp/'
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
        
        comm_world.Barrier()
        
        sim.checkpoint('./validation_cp_mpi/cp/boundiff')

    ########################################################################

    def test_bdiff_n4(self):
        myrank = steps.mpi.rank

        mdl = smodel.Model()
        
        X = smodel.Spec('X', mdl)
        cytosolv = smodel.Volsys('cytosolv', mdl)
        dif_X = smodel.Diff('diffX', cytosolv, X)
        dif_X.setDcst(DCST)
        
        ########################################################################

        mesh = meshio.loadMesh('./validation_rd/meshes/' +MESHFILE)[0]
        
        a = mesh.getBoundMax()[2]-mesh.getBoundMin()[2]
        area = mesh.getMeshVolume()/a
        
        ntets = mesh.countTets()
        comp = stetmesh.TmComp('cyto', mesh, range(ntets))
        comp.addVolsys('cytosolv')
        
        self.assertTrue(SAMPLE == ntets)
        
        numfilled = 0
        while (numfilled < SAMPLE):
            tetidxs[numfilled] = numfilled
            numfilled +=1
        
        # Now find the distance of the center of the tets to the Z lower face
        for i in range(SAMPLE):
            baryc = mesh.getTetBarycenter(int(tetidxs[i]))
            min = mesh.getBoundMin()
            r = baryc[2] - min[2]
            # Convert to microns
            tetrads[i] = r*1.0e6
        
        # And fetch the total number of tets to make the data structures
        ntets = mesh.countTets()

        rng = srng.create('mt19937', 512) 
        rng.initialize(int(time.time()%4294967295)) # The max unsigned long
        
        sim = ssolv.TetVesicle(mdl, mesh, rng, ssolv.EF_NONE)
        sim.setOutputSync(True, 0);
        sim.reset()
        
        sim.setVesicleDT(DT)

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
            tettemp = mesh.getTetTetNeighb(i)
            templist = [t for t in range(4) if tettemp[t] == UNKNOWN_TET]
            if templist:
                boundtets.append(i)
                bt_srftriidx.append(templist)

        self.assertTrue(len(boundtets) == len(bt_srftriidx))

        minztets = []
        boundminz = mesh.getBoundMin()[2] + 0.01e-06
        num2s=0
        for i in range(boundtets.__len__()):
            # get the boundary triangle
            if (bt_srftriidx[i].__len__() == 2): num2s+=1
            for btriidx in bt_srftriidx[i]:
                zminboundtri = True
                tribidx = mesh.getTetTriNeighb(boundtets[i])[btriidx]
                tritemp = mesh.getTri(tribidx)
                trizs = [0.0, 0.0, 0.0]
                trizs[0] = mesh.getVertex(tritemp[0])[2]
                trizs[1] = mesh.getVertex(tritemp[1])[2]
                trizs[2] = mesh.getVertex(tritemp[2])[2]
                for j in range(3):
                    if (trizs[j]>boundminz): zminboundtri = False
                if (zminboundtri): minztets.append(boundtets[i])

        nztets = minztets.__len__()
        volztets = 0.0
        for z in minztets:
            volztets += mesh.getTetVol(z)
        conc = NITER*6.022e23*1.0e-3/volztets

        seed = int(time.time()%4294967295)
        for j in range(NITER):
            sim.reset()
            sim.restore('./validation_cp_mpi/cp/boundiff')
            rng.initialize(seed)
            seed += 1
            for i in range(ntpnts):
                sim.run(tpnts[i])
                for k in range(SAMPLE):
                    res[j, i, k] = sim.getTetSpecCount(int(tetidxs[k]), 'X')

        itermeans = numpy.mean(res, axis = 0)

        ########################################################################

        D = DCST
        pi = math.pi
        nmax = 1000
        N = NINJECT
        N = int((1.0*NINJECT)/nztets)*nztets
        def getprob(x,t):
            if(x>a): 
                print('x out of bounds')
                return
            p=0.0
            for n in range(nmax):
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
            for i in range(NBINS):
                bin_counts[i] = []
            filled = 0
            
            for i in range(itermeans[t].size):
                irad = tetrads[i]
                
                for b in range(NBINS):
                    if(irad>=binmins[b] and irad<binmins[b+1]):
                        bin_counts[b].append(itermeans[t][i])
                        bin_vols[b]+=sim.getTetVol(int(tetidxs[i]))
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
                    self.assertTrue(tol_funcs.tolerable(det_conc, steps_conc, tolerance))


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestBDiffTetvesicleCP, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

