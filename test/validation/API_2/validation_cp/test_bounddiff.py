####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

########################################################################

# 1D diffusion in a finite tube with fully-reflective boundary.
# Circular "point source" from one face. 
# RESTORE

# AIMS: to verify checkpointing and restoring of the spatial stochastic 
# solver 'Tetexact' in the context of the Bounded Diffusion model 
# (see validation/bounddiff.py)
  
########################################################################

import unittest

import steps.interface

from steps.model import *
from steps.sim import *
from steps.geom import *
from steps.rng import *
from steps.saving import *

import os
import math
import time
import numpy

from scipy.constants import Avogadro
from . import tol_funcs
    
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

class TestBDiff(unittest.TestCase):

    def setUp(self):
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

        self.assertTrue(SAMPLE == ntets)
            
        ########################################################################

        rng = RNG('mt19937', 512, int(time.time()%4294967295))
        sim = Simulation('Tetexact', mdl, mesh, rng)

        boundminz = mesh.bbox.min.z + 0.01e-06
        minztets = TetList(mesh=mesh)
        for tri in mesh.surface:
            if all(v.z <= boundminz for v in tri.verts):
                minztets |= tri.tetNeighbs

        sim.newRun()

        sim.TETS(minztets).X.Count = int(NINJECT/len(minztets))

        new_dir = './validation_cp/cp/'
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
            
        sim.checkpoint('./validation_cp/cp/boundiff')

    def test_bdiff(self):

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

        self.assertTrue(SAMPLE == ntets)

        # create the array of tet indices to be found at random
        sampleTets = mesh.tets[:SAMPLE]
        # Now find the distance of the center of the tets to the Z lower face
        minz = mesh.bbox.min.z
        tetrads = [(tet.center.z - minz) * 1e6 for tet in sampleTets]

        minztets = TetList(mesh=mesh)
        for tri in mesh.surface:
            if all(v.z <= minz + 0.01e-6 for v in tri.verts):
                minztets |= tri.tetNeighbs
        nztets = len(minztets)
        
        ########################################################################
        
        rng = RNG('mt19937', 512, int(time.time()%4294967295))
        sim = Simulation('Tetexact', mdl, mesh, rng)

        rs = ResultSelector(sim)

        res = rs.TETS(sampleTets).X.Count

        sim.toSave(res, dt=DT)

        new_dir = './validation_cp/cp/'
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)

        for j in range(NITER):
            sim.newRun()
            sim.restore('./validation_cp/cp/boundiff')
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

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestBDiff, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

