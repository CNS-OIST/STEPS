########################################################################

# Stochastic second-order irreversible reaction on a surface.
# RESTORE

# AIMS: to verify checkpointing and restoring of the well-mixed stochastic 
# solver 'Wmdirect' in the context of the Second-Order Irreversible 
# Surface Reaction model (see validation/second_order_irev_2D.py)
  
########################################################################

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import math
import time 
import numpy

from tol_funcs import *

print "Reaction - Second-order, irreversible, 2D:"
import second_order_irev_2D_cp

########################################################################

def test_soirev2d():
    VOL = 1.0e-18

    COUNTA = 100.0
    n=2.0
    COUNTB = COUNTA/n 


    KCST = 10.0e10			# The reaction constant

    AREA = 10.0e-12

    CCST = KCST/(6.02214179e23*AREA)


    NITER = 1000			# The number of iterations
    DT = 0.05			# Sampling time-step
    INT = 1.05			# Sim endtime

    # In tests fewer than 0.1% fail with tolerance of 2%
    tolerance = 2.0/100

    ########################################################################

    mdl  = smod.Model()

    A = smod.Spec('A', mdl)
    B = smod.Spec('B', mdl)
    C = smod.Spec('C', mdl)

    surfsys = smod.Surfsys('ssys',mdl)

    SR1 = smod.SReac('SR1', surfsys, slhs = [A, B], srhs = [C], kcst = KCST)

    geom = sgeom.Geom()

    comp1 = sgeom.Comp('comp1', geom, VOL)
    patch1 = sgeom.Patch('patch1', geom, comp1, area = AREA)
    patch1.addSurfsys('ssys')

    import random
    rng = srng.create('mt19937', 1000)
    rng.initialize(int(random.random()*4294967295))


    sim = ssolv.Wmdirect(mdl, geom, rng)
    sim.reset()

    tpnts = numpy.arange(0.0, INT, DT)
    ntpnts = tpnts.shape[0]

    res_m = numpy.zeros([NITER, ntpnts, 3])

    for i in range (0, NITER):
        sim.restore('./validation_cp/cp/second_order_irev_2D')
        for t in range(0, ntpnts):
            sim.run(tpnts[t])
            res_m[i, t, 0] = sim.getPatchCount('patch1', 'A')
            res_m[i, t, 1] = sim.getPatchCount('patch1', 'B')       

    mean_res = numpy.mean(res_m, 0)

    lnBA = numpy.zeros(ntpnts)
    lineAB = numpy.zeros(ntpnts)

    C = COUNTA-COUNTB
    passed  =True
    max_err = 0.0

    for i in range(ntpnts):
        A = mean_res[i][0]
        B = mean_res[i][1]
        lnBA[i] = math.log(B/A)
        lineAB[i] = math.log(COUNTB/COUNTA) -C*CCST*tpnts[i]
        assert(tolerable(lnBA[i], lineAB[i], tolerance))


