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

import unittest

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import math
import time 
import numpy
import steps.utilities.meshio as meshio

from . import tol_funcs

########################################################################
class TestRDKisilevich(unittest.TestCase):

    def test_kisilevich(self):
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

        # create the array of tet indices to be found at random
        tetidxs = numpy.zeros(SAMPLE, dtype = 'int')
        # further create the array of tet barycenter distance to center
        tetrads = numpy.zeros(SAMPLE)

        mdl  = smod.Model()

        A = smod.Spec('A', mdl)
        B = smod.Spec('B', mdl)

        volsys = smod.Volsys('vsys',mdl)

        R1 = smod.Reac('R1', volsys, lhs = [A,B], rhs = [])

        R1.setKcst(RCST)

        D_a = smod.Diff('D_a', volsys, A)
        D_a.setDcst(DCSTA)
        D_b = smod.Diff('D_b', volsys, B)
        D_b.setDcst(DCSTB)

        mesh = meshio.loadMesh('validation_rd/meshes/brick_40_4_4_1686tets')[0]

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
        self.assertTrue(SAMPLE <= ntets)

        numfilled = 0
        while (numfilled < SAMPLE):
            tetidxs[numfilled] = numfilled
            numfilled +=1

        # Now find the distance of the center of the tets to the Z lower face
        for i in range(SAMPLE):
            baryc = mesh.getTetBarycenter(int(tetidxs[i]))
            r = baryc[0]
            tetrads[i] = r*1.0e6

        Atets = acomptets
        Btets = bcomptets

        rng = srng.create('r123', 512)
        rng.initialize(1000)

        sim = ssolv.Tetexact(mdl, mesh, rng)

        sim.reset()

        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]

        resA = numpy.zeros((NITER, ntpnts, SAMPLE))
        resB = numpy.zeros((NITER, ntpnts, SAMPLE))


        for i in range (0, NITER):
            sim.reset()
            
            sim.setDiffBoundarySpecDiffusionActive('diffb', 'A', True)
            sim.setDiffBoundarySpecDiffusionActive('diffb', 'B', True)
            
            sim.setCompSpecCount('compa', 'A', NA0)
            sim.setCompSpecCount('compb', 'B', NB0)
            
            for t in range(0, ntpnts):
                sim.run(tpnts[t])
                for k in range(SAMPLE):
                    resA[i,t,k] = sim.getTetSpecCount(int(tetidxs[k]), 'A')
                    resB[i,t,k] = sim.getTetSpecCount(int(tetidxs[k]), 'B')


        itermeansA = numpy.mean(resA, axis=0)
        itermeansB = numpy.mean(resB, axis=0)

        def getdetc(t, x):
            N = 1000        # The number to represent infinity in the exponential calculation
            L = 20e-6
            
            concA  = 0.0
            for n in range(N):
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
            
            for b in range(NBINS+1):
                binmins[b] = r
                if (b!=NBINS): tetradsbinned[b] = r +rsec/2.0
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
                        bin_vols[b]+=sim.getTetVol(int(tetidxs[i]))
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
                bin_concsA[c]*=(1.0e-3/6.022e23)*1.0e6    
                bin_concsB[c]/=(bin_vols[c])
                bin_concsB[c]*=(1.0e-3/6.022e23)*1.0e6 
            
            for i in range(NBINS):
                rad = abs(tetradsbinned[i])*1.0e-6
                
                if (tetradsbinned[i] < -5):
                    # compare A
                    det_conc = getdetc(tpnts[tidx], rad)
                    steps_conc = bin_concsA[i]
                    self.assertTrue(tol_funcs.tolerable(det_conc, steps_conc, tolerance))

                if (tetradsbinned[i] > 5):
                    # compare B
                    det_conc = getdetc(tpnts[tidx], rad)
                    steps_conc = bin_concsB[i]
                    self.assertTrue(tol_funcs.tolerable(det_conc, steps_conc, tolerance))

    ########################################################################

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestRDKisilevich))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

