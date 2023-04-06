########################################################################

# 2D diffusion on an infinite strip from line "point source". 

# AIMS: to verify STEPS spatial-stochastic solver 'Tetexact' supports 
# local initial conditions and calculates surface diffusion rates correctly.

# STEPS equivalent model: Stochastic 2D diffusion on outer surface of 
# ring; source in small region of x=10um triangles.

# Verification also takes place of model and mesh construction 
# components, particularly mesh loading and manipulation capabilities 
# with functions such as steps.geom.Tetmesh.getTriBarycenter and 
# steps.geom.Tetmesh.getTriArea etc. 
# Localised recording by steps.solver.Tetexact.getTriCount is also verified. 

# A 5% tolerance is imposed when comparing the mean output from 100 
# stochastic simulations of the STEPS model to the analytical solution. 
# There is an expected probability of failure of < 1%.
  
########################################################################

import unittest

import steps.interface

from steps.model import *
from steps.sim import *
from steps.saving import *
from steps.geom import *
from steps.rng import *

import datetime
import time
import math
import numpy

from . import tol_funcs

########################################################################

class TestRDUnbdiff2DLinesourceRing(unittest.TestCase):

    def setUp(self):
        global NITER, DT, INT, NINJECT, DCST, tolerance, MESHFILE

        # Number of iterations; plotting dt; sim endtime:
        NITER = 100
        DT = 0.02
        INT = 3.00

        # Number of molecules injected in center; diff constant
        NINJECT = 5000
        DCST = 0.01e-9

        # In tests <1% fail with tolerance of 5%
        tolerance = 5.0/100

        ########################################################################

        MESHFILE = 'ring2_or10_ir9_injx0width2_640tets.inp'

    ########################################################################


    ########################################################################

    def test_unbdiff2D_linesource_ring(self):
        "Surface Diffusion - Unbounded, line source (Tetexact)"

        mdl = Model()
        with mdl:
            X = Species.Create()
            ssys = SurfaceSystem.Create()
            with ssys:
                dif_X = Diffusion.Create(X, DCST)
        
        ########################################################################

        mesh = TetMesh.Load('validation_rd/meshes/'+MESHFILE)
            
        ntets = len(mesh.tets)
        with mesh:
            comp = Compartment.Create(mesh.tets)

            patch_tris = TriList(tri for tri in mesh.surface if 0.00000995 < numpy.linalg.norm(tri.center[0:2]) < 0.00001001)

            patch = Patch.Create(patch_tris, comp, None, ssys)
        
            # X should be maximum (10) for the inject region
            inject_tris = TriList(tri for tri in patch_tris if tri.center.x > 9.999e-6)

        sample_tris = patch_tris - inject_tris
                
        # Now find the distances along the edge for all tris
        tridists = []
        triareas = []
        for tri in sample_tris:
            # Triangles will be separated into those to the 'high' side (positive y) with positive L
            # and those on the low side (negative y) with negative L 
            rad = numpy.linalg.norm(tri.center[0:2])
            theta = math.atan2(tri.center.y, tri.center.x)
            tridists.append(theta * rad* 1e6)
            triareas.append(tri.Area)

        rng = RNG('r123', 1024, 1000)

        sim = Simulation('Tetexact', mdl, mesh, rng)

        rs = ResultSelector(sim)

        res_count = rs.TRIS(sample_tris).X.Count
        res_conc = 1e-12 * rs.TRIS(sample_tris).X.Count / rs.TRIS(sample_tris).Area

        sim.toSave(res_count, res_conc, dt=DT)

        for j in range(NITER):
            sim.newRun()

            sim.TRIS(inject_tris).X.Count = float(NINJECT) / len(inject_tris)

            sim.run(INT)

        itermeans_count = numpy.mean(res_count.data, axis = 0)
        itermeans_conc = numpy.mean(res_conc.data, axis = 0)

        ########################################################################

        tpnt_compare = [100, 150]

        passed = True
        max_err = 0.0

        for t in tpnt_compare:
            bin_n = 50
            
            r_min=0
            r_max=0
            
            for i in tridists: 		
                if (i > r_max):
                    r_max = i
                if (i < r_min):
                    r_min = i
            
            r_seg = (r_max-r_min)/bin_n
            bin_mins = numpy.zeros(bin_n+1)
            r_tris_binned = numpy.zeros(bin_n)
            bin_areas = numpy.zeros(bin_n)    
            
            r = r_min
            for b in range(bin_n + 1):
                bin_mins[b] = r
                if (b!=bin_n):
                    r_tris_binned[b] = r +r_seg/2.0
                r+=r_seg
            bin_counts = [None]*bin_n
            for i in range(bin_n): bin_counts[i] = []
            for i in range((itermeans_count[t].size)):
                i_r = tridists[i]
                for b in range(bin_n):
                    if(i_r>=bin_mins[b] and i_r<bin_mins[b+1]):
                        bin_counts[b].append(itermeans_count[t][i])
                        bin_areas[b]+=triareas[i]
                        break
            
            bin_concs = numpy.zeros(bin_n)
            for c in range(bin_n): 
                for d in range(bin_counts[c].__len__()):
                    bin_concs[c] += bin_counts[c][d]
                bin_concs[c]/=(bin_areas[c]*1.0e12)
            
            for i in range(bin_n):
                if (r_tris_binned[i] > -10.0 and r_tris_binned[i] < -2.0) \
                or (r_tris_binned[i] > 2.0 and r_tris_binned[i] < 10.0):
                    dist = r_tris_binned[i]*1e-6
                    det_conc = 1e-6*(NINJECT/(4*math.sqrt((math.pi*DCST*res_count.time[0,t]))))*(math.exp((-1.0*(dist*dist))/(4*DCST*res_count.time[0,t])))	
                    steps_conc = bin_concs[i]
                    self.assertTrue(tol_funcs.tolerable(det_conc, steps_conc, tolerance))

    ########################################################################

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestRDUnbdiff2DLinesourceRing, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
