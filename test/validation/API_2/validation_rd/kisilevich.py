########################################################################

# Stochastic degradation-diffusion process.

# AIMS: to verify STEPS spatial stochastic solver 'Tetexact' 
# supports spatially-separated initial conditions and computes
# diffusion and reaction rates correctly when applied to a 
# degradation-diffusion process. Also verifies performance of the 
# "Diffusion boundary" between two different spatial compartments.

# For a more detailed description of the analytical system and 
# equivalent STEPS model see:
# http://www.biomedcentral.com/content/supplementary/1752-0509-6-36-s4.pdf
# "Degradation-diffusion process with initially separated reactants"

# Verification also takes place of the necessary steps to build the model, 
# such as multiple mesh-based compartment creation, diffusion boundary 
# activation, and recording from tetrahedrons in different compartments. 

# A 7.5% tolerance is imposed when comparing the mean output from 50 
# stochastic simulations of the STEPS model to the analytical solution. 
# There is an expected probability of failure of < 1%.
  
########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import math
import time 
import numpy

from scipy.constants import Avogadro
from . import tol_funcs

########################################################################

def test_kisilevich():
    "Reaction-diffusion - Degradation-diffusion (Tetexact)"

    NITER = 50       # The number of iterations
    DT = 0.1         # Sampling time-step
    INT = 0.3        # Sim endtime

    DCSTA = 400*1e-12
    DCSTB = DCSTA
    RCST = 100000.0e6

    #NA0 = 100000    # 1000000            # Initial number of A molecules
    NA0 = 1000
    NB0 = NA0        # Initial number of B molecules

    SAMPLE = 1686

    # <1% fail with a tolerance of 7.5%
    tolerance = 7.5/100

    ########################################################################
    mdl = Model()
    r = ReactionManager()
    with mdl:
        SA, SB = Species.Create()
        volsys = VolumeSystem.Create()
        with volsys:
            SA + SB >r['R1']> None
            r['R1'].K = RCST
            
            D_a =     Diffusion.Create(SA, DCSTA)
            D_b =     Diffusion.Create(SB, DCSTB)


    mesh = TetMesh.Load('./validation_rd/meshes/brick_40_4_4_1686tets')

    with mesh:
        acomptets = TetList(tet for tet in mesh.tets if tet.center.x < 0)
        bcomptets = mesh.tets - acomptets
        
        compa = Compartment.Create(acomptets, volsys)
        compb = Compartment.Create(bcomptets, volsys)
        
        diffb = DiffBoundary.Create(acomptets.surface & bcomptets.surface)

    # create the array of tet indices to be found at random
    sampleTets = mesh.tets[:SAMPLE]
    # Now find the distance of the center of the tets to the Z lower face
    tetrads = [tet.center.x * 1e6 for tet in sampleTets]

    VOLA = mesh.Vol/2.0
    VOLB = VOLA

    rng = RNG('r123', 512, 1000)

    sim = Simulation('Tetexact', mdl, mesh, rng)

    rs = ResultSelector(sim)

    resA = rs.TETS(sampleTets).SA.Count
    resB = rs.TETS(sampleTets).SB.Count

    sim.toSave(resA, resB, dt=DT)

    for i in range (0, NITER):    
        sim.newRun()

        sim.diffb.SA.DiffusionActive = True
        sim.diffb.SB.DiffusionActive = True
        sim.compa.SA.Count = NA0
        sim.compb.SB.Count = NB0

        sim.run(INT)

    itermeansA = numpy.mean(resA.data, axis=0)
    itermeansB = numpy.mean(resB.data, axis=0)

    def getdetc(t, x):
        N = 1000        # The number to represent infinity in the exponential calculation
        L = 20e-6
        
        concA  = 0.0
        for n in range(N):
            concA+= ((1.0/(2*n +1))*math.exp((-(DCSTA/(20.0e-6))*math.pow((2*n +1), 2)*math.pow(math.pi, 2)*t)/(4*L))*math.sin(((2*n +1)*math.pi*x)/(2*L)))
        concA*=((4*NA0/math.pi)/(VOLA*Avogadro * 1e3))*1.0e6    
        
        return concA


    tpnt_compare = [1, 2]
    passed = True
    max_err = 0.0

    for tidx in tpnt_compare:
        NBINS=10
        radmax = max(tetrads)
        radmin = min(tetrads)
        
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
        
        bin_countsA = [None]*NBINS
        bin_countsB = [None]*NBINS
        for i in range(NBINS):
            bin_countsA[i] = []
            bin_countsB[i] = []
        filled = 0
        
        for i in range(itermeansA[tidx].size):
            irad = tetrads[i]
            
            for b in range(NBINS):
                if(irad>=binmins[b] and irad<binmins[b+1]):
                    bin_countsA[b].append(itermeansA[tidx][i])
                    bin_vols[b]+=sampleTets[i].Vol
                    filled+=1.0
                    break
        filled = 0
        for i in range(itermeansB[tidx].size):
            irad = tetrads[i]
            
            for b in range(NBINS):
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
            bin_concsA[c]*=(1.0e-3/Avogadro)*1.0e6    
            bin_concsB[c]/=(bin_vols[c])
            bin_concsB[c]*=(1.0e-3/Avogadro)*1.0e6 
        
        for i in range(NBINS):
            rad = abs(tetradsbinned[i])*1.0e-6
            
            if (tetradsbinned[i] < -5):
                # compare A
                det_conc = getdetc(resA.time[0,tidx], rad)
                steps_conc = bin_concsA[i]
                assert tol_funcs.tolerable(det_conc, steps_conc, tolerance)

            if (tetradsbinned[i] > 5):
                # compare B
                det_conc = getdetc(resB.time[0,tidx], rad)
                steps_conc = bin_concsB[i]
                assert tol_funcs.tolerable(det_conc, steps_conc, tolerance)

########################################################################
# END

