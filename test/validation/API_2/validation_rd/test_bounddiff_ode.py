########################################################################

# 1D diffusion in a finite tube with fully-reflective boundary.
# Circular "point source" from one face. 

# AIMS: to verify STEPS spatial-deterministic solver 'TetODE' supports 
# local initial condition on boundary, calculates volume diffusion 
# rates correctly and imposes reflective boundary condition.

# STEPS equivalent model: Deterministic 3D diffusion in a cylinder; source 
# from tetrahedrons on boundary of one circular face.

# For a more detailed description of the analytical system and 
# equivalent STEPS model see:
# http://www.biomedcentral.com/content/supplementary/1752-0509-6-36-s4.pdf
# "1D diffusion in a finite tube from a point source at one end"

# Verification also takes place of model and mesh construction 
# components, particularly mesh loading and manipulation capabilities 
# with functions such as steps.utilities.meshio.loadMesh and 
# steps.geom.Tetmesh.getTetTriNeighb, steps.geom.Tetmesh.getTri etc. 
# Localised recording by steps.solver.TetODE.getTetCount is also verified. 

# Even though this is a deterministic model, a tolerance of 2.5% is 
# permitted. This is because an infinitely thin plane source is 
# not replicated in STEPS and a small error is introduced 
# by the initial variance introduced by the tetrahedral source.
  
########################################################################

import unittest

import steps.interface

from steps.geom import *
from steps.model import *
from steps.sim import *
from steps.saving import *
from steps.geom import *
from steps.rng import *

import math
import time
import numpy

from scipy.constants import Avogadro
from . import tol_funcs

########################################################################

class TestRDBoundDiffODE(unittest.TestCase):

    def setUp(self):
        global NITER, DT, INT, NINJECT, DCST, tolerance, SAMPLE, MESHFILE, tetidxs, tetrads

        NITER = 1
        DT = 0.01
        INT = 0.10

        # The number of initial molecules:
        NINJECT = 10000	

        DCST = 0.2e-9

        # Small expected error from non plane source
        tolerance = 2.5/100

        # The number of tets to sample at random:
        SAMPLE = 1060	

        MESHFILE = 'cyl_diam2__len10_1060tets'

        # create the array of tet indices to be found at random
        tetidxs = numpy.zeros(SAMPLE, dtype = 'int')
        # further create the array of tet barycenter distance to center
        tetrads = numpy.zeros(SAMPLE)



    ########################################################################

    def test_bounddiff_ode(self):
        "Diffusion - Bounded (TetODE)"

        ########################################################################
        mdl = Model()
        with mdl:
            X = Species.Create()
            cytosolv = VolumeSystem.Create()
            with cytosolv:
                dif_X = Diffusion.Create(X, DCST)

        ########################################################################

        mesh = TetMesh.Load('validation_rd/meshes/' + MESHFILE)
        
        a = mesh.bbox.max.z - mesh.bbox.min.z
        area = mesh.Vol/a
        
        ntets = len(mesh.tets)
        with mesh:
            comp = Compartment.Create(range(ntets), 'cytosolv')
        
        self.assertTrue(SAMPLE == ntets)

        # create the array of tet indices to be found at random
        sampleTets = mesh.tets[:SAMPLE]
        # Now find the distance of the center of the tets to the Z lower face
        minz = mesh.bbox.min.z
        tetrads = [(tet.center.z - minz) * 1e6 for tet in sampleTets]

        # Find the tets connected to the bottom face
        minztets = TetList(mesh=mesh)
        for tri in mesh.surface:
            if all(v.z <= minz + 0.01e-6 for v in tri.verts):
                minztets |= tri.tetNeighbs
        nztets = len(minztets)

        sim = Simulation('TetODE', mdl, mesh, None)
        sim.setTolerances(1e-3, 1e-3)

        rs = ResultSelector(sim)

        res = rs.TETS(sampleTets).X.Count

        sim.toSave(res, dt=DT)

        for j in range(NITER):
            sim.newRun()
            sim.TETS(minztets).X.Count = int((1.0*NINJECT)/nztets)
            sim.run(INT)

        itermeans = numpy.mean(res.data, axis = 0)

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
                if (n==0):
                    A = math.sqrt(1.0/a)
                else :
                    A = math.sqrt(2.0/a)
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
                if (r > radmax):
                    radmax = r
                if (r < radmin) :
                    radmin = r
            
            rsec = (radmax-radmin)/NBINS
            binmins = numpy.zeros(NBINS+1)
            tetradsbinned = numpy.zeros(NBINS)
            r = radmin
            bin_vols = numpy.zeros(NBINS)
            
            for b in range(NBINS+1):
                binmins[b] = r
                if (b!=NBINS):
                    tetradsbinned[b] = r +rsec/2.0
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
                        bin_vols[b]+=sampleTets[i].Vol
                        filled+=1.0
                        break
            bin_concs = numpy.zeros(NBINS)
            for c in range(NBINS): 
                for d in range(bin_counts[c].__len__()):
                    bin_concs[c] += bin_counts[c][d]
                bin_concs[c]/=(bin_vols[c])
                bin_concs[c]*=(1.0e-3/Avogadro)*1.0e6
            
            for i in range(NBINS):
                if (tetradsbinned[i] > 2 and tetradsbinned[i] < 8):
                    rad = tetradsbinned[i]*1.0e-6
                    det_conc = (getprob(rad, res.time[0,t])/area)*(1.0/(Avogadro * 1e-3))
                    steps_conc = bin_concs[i]
                    self.assertTrue(tol_funcs.tolerable(det_conc, steps_conc, tolerance))

    ########################################################################

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestRDBoundDiffODE, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
