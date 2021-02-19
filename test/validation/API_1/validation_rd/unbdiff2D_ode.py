########################################################################

# 2D diffusion on an infinite surface from point source. 

# AIMS: to verify STEPS spatial-deterministic solver 'TetODE' supports 
# local initial conditions and calculates surface diffusion rates correctly.

# STEPS equivalent model: Deterministic 2D diffusion on one outer surface 
# of a disk; central triangle source.

# This model is equivalet to that described at
# http://www.biomedcentral.com/content/supplementary/1752-0509-6-36-s4.pdf
# "3D diffusion in infinite volume"
# with a 2D compartments ("patch") instead of a 3D compartment.

# Verification also takes place of model and mesh construction 
# components, particularly mesh loading and manipulation capabilities 
# with functions such as steps.geom.Tetmesh.getSurfTris and 
# steps.geom.Tetmesh.getTriArea etc. 
# Localised recording by steps.solver.TetODE.getTriCount is also verified. 

# Even though this is a deterministic model, a tolerance of 2.5% is 
# permitted. This is because a point source is not replicated in STEPS 
# and a small error is introduced by the initial variance from the  
# finite triangle source.

########################################################################

import steps.model as smodel
import steps.solver as solvmod
import steps.utilities.meshio as smeshio
import steps.geom as stetmesh
import steps.rng as srng

import datetime
import time
import math
import numpy

from . import tol_funcs

########################################################################

def setup_module():
    global DT, INT, NINJECT, DCST, tolerance, MESHFILE

    # plotting dt; sim endtime:
    DT = 0.001
    INT = 0.151

    # Number of molecules injected in center; diff constant
    NINJECT = 1000 	

    DCST = 0.02e-9

    # Due to non-point source there is a small expected error, maximum 2.5% 
    tolerance = 2.5/100

    ########################################################################

    MESHFILE = 'coin_10r_1h_99k.inp'

########################################################################

def gen_model():
    mdl = smodel.Model()
    X = smodel.Spec('X', mdl)
    
    ssys = smodel.Surfsys('ssys', mdl)
    dif_X = smodel.Diff('diffX', ssys, X)
    dif_X.setDcst(DCST)
    
    return mdl

########################################################################

def gen_geom():
    mesh = smeshio.loadMesh('validation_rd/meshes/'+MESHFILE)[0]
    
    ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.5e-6])
    
    ntets = mesh.countTets()
    comp = stetmesh.TmComp('cyto', mesh, range(ntets))
    
    alltris = mesh.getSurfTris()
    
    patch_tris = []
    for t in alltris:
        vert0, vert1, vert2 = mesh.getTri(t)
        if (mesh.getVertex(vert0)[2] > 0.0 and mesh.getVertex(vert1)[2] > 0.0 and mesh.getVertex(vert2)[2] > 0.0):
            patch_tris.append(t)
    
    patch_tris_n = len(patch_tris)
                
    patch = stetmesh.TmPatch('patch', mesh, patch_tris, icomp = comp)
    patch.addSurfsys('ssys')
    
    trirads = numpy.zeros(patch_tris_n)
    triareas = numpy.zeros(patch_tris_n)
                
    # Find the central triangle
    ctet_trineighbs = mesh.getTetTriNeighb(ctetidx)
    ctri_idx=-1
    for t in ctet_trineighbs: 
        if t in patch_tris:
            ctri_idx = t
    
    # Now find the distance of the center of the tets to the center of the center tet (at 0,0,0)
    cbaryc = mesh.getTriBarycenter(ctri_idx)
    for i in range(patch_tris_n):
        baryc = mesh.getTriBarycenter(patch_tris[i])
        r2 = math.pow((baryc[0]-cbaryc[0]),2) + math.pow((baryc[1]-cbaryc[1]),2) + math.pow((baryc[2]-cbaryc[2]),2)
        r = math.sqrt(r2)
        # Conver to microns
        trirads[i] = r*1.0e6
        triareas[i] = mesh.getTriArea(patch_tris[i])
        
    return mesh, patch_tris, patch_tris_n, ctri_idx, trirads, triareas

########################################################################

def test_unbdiff2D_ode():
    "Surface Diffusion - Unbounded, point source (TetODE)"

    m = gen_model()
    g, patch_tris, patch_tris_n, ctri_idx, trirads, triareas = gen_geom()

    sim = solvmod.TetODE(m, g)
    sim.setTolerances(1.0e-7, 1.0e-7)

    tpnts = numpy.arange(0.0, INT, DT)
    ntpnts = tpnts.shape[0]

    res_count = numpy.zeros((ntpnts, patch_tris_n))

    sim.setTriCount(ctri_idx, 'X', NINJECT)

    for i in range(ntpnts):
        sim.run(tpnts[i])
        for k in range(patch_tris_n):
            res_count[i, k] = sim.getTriCount(patch_tris[k], 'X')

    ########################################################################

    tpnt_compare = [50, 100, 150]
    passed = True
    max_err = 0.0

    for t in tpnt_compare:
        bin_n = 20
        
        r_max = 0.0 	
        for i in trirads: 		
            if (i > r_max): r_max = i 	
        
        r_min = 0.0
        
        r_seg = (r_max-r_min)/bin_n
        bin_mins = numpy.zeros(bin_n+1)
        r_tris_binned = numpy.zeros(bin_n)
        bin_areas = numpy.zeros(bin_n)    
        
        r = r_min
        for b in range(bin_n + 1):
            bin_mins[b] = r
            if (b!=bin_n): r_tris_binned[b] = r +r_seg/2.0
            r+=r_seg
        bin_counts = [None]*bin_n
        for i in range(bin_n): bin_counts[i] = []
        for i in range((res_count[t].size)):
            i_r = trirads[i]
            for b in range(bin_n):
                if(i_r>=bin_mins[b] and i_r<bin_mins[b+1]):
                    bin_counts[b].append(res_count[t][i])
                    bin_areas[b]+=sim.getTriArea(int(patch_tris[i]))
                    break
        
        bin_concs = numpy.zeros(bin_n)
        for c in range(bin_n): 
            for d in range(bin_counts[c].__len__()):
                bin_concs[c] += bin_counts[c][d]
            bin_concs[c]/=(bin_areas[c]*1.0e12)
        
        for i in range(bin_n):
            if (r_tris_binned[i] > 1.0 and r_tris_binned[i] < 5.0):
                rad = r_tris_binned[i]*1.0e-6
                det_conc = 1.0e-12*(NINJECT/(4*math.pi*DCST*tpnts[t]))*(math.exp((-1.0*(rad*rad))/(4*DCST*tpnts[t])))
                steps_conc = bin_concs[i]
                assert tol_funcs.tolerable(det_conc, steps_conc, tolerance)

########################################################################
# END
