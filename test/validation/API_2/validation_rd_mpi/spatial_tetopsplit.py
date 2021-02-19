########################################################################

# Well-mixed chemical reactions. 

# AIMS: to verify STEPS well-mixed stochastic solver 'TetOpSplit' computes 
# reaction rates correctly when many different chemical species and 
# reactions are present in the model. 

# This model combines all well-mixed validations, as described at: 
# http://www.biomedcentral.com/content/supplementary/1752-0509-6-36-s4.pdf
# and in the equivalent individual model scripts.

########################################################################
 
import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import numpy
import math
import time

from scipy.constants import Avogadro
from . import tol_funcs

########################################################################

def setup_module():
    global ran_sim

    ran_sim = False

def run_sim():
    # Set up and run the simulations once, before the tests
    # analyze the results.

    ##################### First order irreversible #########################

    global KCST_foi, N_foi, tolerance_foi

    KCST_foi = 5              # The reaction constant
    N_foi = 50                # Can set count or conc

    NITER_foi = 100000        # The number of iterations

    # Tolerance for the comparison:
    # In test runs, with good code, < 1%  will fail with a 1.5% tolerance
    tolerance_foi = 2.0/100  

    ####################### First order reversible #########################

    global KCST_f_for, KCST_b_for, COUNT_for, tolerance_for

    KCST_f_for = 10.0         # The reaction constant
    KCST_b_for = 2.0

    COUNT_for = 100000        # Can set count or conc

    NITER_for = 10            # The number of iterations

    # In test runs, with good code, <0.1% will fail with a tolerance of 1% 
    tolerance_for = 1.0/100

    ####################### Second order irreversible A2 ###################

    global KCST_soA2, CONCA_soA2, tolerance_soA2

    KCST_soA2 = 10.0e6        # The reaction constant

    CONCA_soA2 = 10.0e-6

    NITER_soA2 = 1000         # The number of iterations

    # In test runs, with good code, <0.1% will fail with a tolerance of 2%
    tolerance_soA2 = 3.0/100

    ####################### Second order irreversible AA ###################

    global KCST_soAA, CONCA_soAA, CONCB_soAA, tolerance_soAA
    
    KCST_soAA = 5.0e6        # The reaction constant
    
    CONCA_soAA = 20.0e-6
    CONCB_soAA = CONCA_soAA
    
    NITER_soAA = 1000         # The number of iterations
    
    # In test runs, with good code, <0.1% will fail with a tolerance of 1%
    tolerance_soAA = 2.0/100
    
    ####################### Second order irreversible AB ###################

    global KCST_soAB, CONCA_soAB, CONCB_soAB, tolerance_soAB

    KCST_soAB = 5.0e6         # The reaction constant

    CONCA_soAB = 1.0e-6
    n_soAB = 2
    CONCB_soAB = CONCA_soAB/n_soAB

    NITER_soAB = 1000         # The number of iterations

    # In test runs, with good code, <0.1% will fail with a tolerance of 1% 
    tolerance_soAB = 1.0/100

    ####################### Third order irreversible A3 ###################

    global KCST_toA3, CONCA_toA3, tolerance_toA3
    
    KCST_toA3 = 1.0e12        # The reaction constant
    
    CONCA_toA3 = 10.0e-6
    
    NITER_toA3 = 1000         # The number of iterations
    
    # In test runs, with good code, <0.1% will fail with a tolerance of 1%
    tolerance_toA3 = 3.0/100

    ####################### Third order irreversible A2B ###################

    global KCST_toA2B, CONCA_toA2B, CONCB_toA2B, tolerance_toA2B
    
    KCST_toA2B = 0.1e12        # The reaction constant
    
    CONCA_toA2B = 30.0e-6
    CONCB_toA2B = 20.0e-6
    
    NITER_toA2B = 1000         # The number of iterations
    
    # In test runs, with good code, <0.1% will fail with a tolerance of 1%
    tolerance_toA2B = 1.0/100
        
    ####################### Second order irreversible 2D ###################

    global COUNTA_so2d, COUNTB_so2d, CCST_so2d, tolerance_so2d

    COUNTA_so2d = 100.0
    n_so2d=2.0
    COUNTB_so2d = COUNTA_so2d/n_so2d 

    KCST_so2d = 10.0e10       # The reaction constant

    AREA_so2d = 10.0e-12


    NITER_so2d = 1000         # The number of iterations

    # In tests fewer than 0.1% fail with tolerance of 2%
    tolerance_so2d = 2.0/100

    ############################ Common parameters ########################

    global VOL

    DT = 0.1                  # Sampling time-step
    INT = 1.1                 # Sim endtime

    NITER_max = 100000

    ########################################################################

    mdl = Model()
    r = ReactionManager()
    with mdl:
        A_foi, A_for, B_for, A_soA2, C_soA2, A_soAA, B_soAA, C_soAA, A_soAB, B_soAB, \
        C_soAB, A_toA3, C_toA3, A_toA2B, B_toA2B, C_toA2B, A_so2d, B_so2d, C_so2d = Species.Create()

        volsys = VolumeSystem.Create()
        surfsys = SurfaceSystem.Create()

        with volsys:
            # First order irreversible
            A_foi_diff = Diffusion.Create(A_foi, 1e-14)
            A_foi >r['R1_foi']> None
            r['R1_foi'].K = KCST_foi

            # First order reversible
            A_for_diff = Diffusion.Create(A_for, 1e-14)
            B_for_diff = Diffusion.Create(B_for, 1e-14)
            A_for <r['R1_for']> B_for
            r['R1_for'].K = KCST_f_for, KCST_b_for

            # Second order irreversible A2
            A_soA2_diff = Diffusion.Create(A_soA2, 1e-12)
            2 * A_soA2 >r['R1_soA2']> C_soA2
            r['R1_soA2'].K = KCST_soA2

            # Second order irreversible AA
            A_soAA_diff = Diffusion.Create(A_soAA, 2e-13)
            B_soAA_diff = Diffusion.Create(B_soAA, 2e-13)
            A_soAA + B_soAA >r['R1_soAA']> C_soAA
            r['R1_soAA'].K = KCST_soAA

            # Second order irreversible AB
            A_soAB_diff = Diffusion.Create(A_soAB, 1e-13)
            B_soAB_diff = Diffusion.Create(B_soAB, 1e-13)
            A_soAB + B_soAB >r['R1_soAB']> C_soAB
            r['R1_soAB'].K = KCST_soAB

            # Third order irreversible A3
            A_soA3_diff = Diffusion.Create(A_toA3, 2e-13)
            3 * A_toA3 >r['R1_toA3']> C_toA3
            r['R1_toA3'].K = KCST_toA3

            # Third order irreversible A2B
            A_soA2B_diff = Diffusion.Create(A_toA2B, 1e-13)
            B_soA2B_diff = Diffusion.Create(B_toA2B, 1e-13)
            2 * A_toA2B + B_toA2B >r['R1_toA2B']> C_toA2B
            r['R1_toA2B'].K = KCST_toA2B
                
        with surfsys:
            # Second order irreversible 2D
            A_so2d_diff = Diffusion.Create(A_so2d, 1e-12)
            B_so2d_diff = Diffusion.Create(B_so2d, 1e-12)
            A_so2d.s + B_so2d.s >r['SR1_so2d']> C_so2d.s
            r['SR1_so2d'].K = KCST_so2d

    mesh = TetMesh.LoadAbaqus('validation_rd_mpi/meshes/sphere_rad1_37tets.inp', scale=1e-06)
    VOL = mesh.Vol
    with mesh:
        comp1 = Compartment.Create(mesh.tets, volsys)
        patch1 = Patch.Create(mesh.surface, comp1, None, surfsys)
    
    CCST_so2d = KCST_so2d/(Avogadro*patch1.Area)

    rng = RNG('r123', 512, 100)

    part = LinearMeshPartition(mesh, MPI.nhosts, 1, 1)

    sim = Simulation('TetOpSplit', mdl, mesh, rng, MPI.EF_NONE, part)

    rs = ResultSelector(sim)

    res_m_foi = rs.comp1.A_foi.Count

    res_m_for = rs.comp1.LIST(A_for, B_for).Conc * 1e6

    res_m_soA2 = rs.comp1.A_soA2.Conc

    res_m_soAA = rs.comp1.LIST(A_soAA, B_soAA).Conc

    res_m_soAB = rs.comp1.LIST(A_soAB, B_soAB).Conc

    res_m_toA3 = rs.comp1.A_toA3.Conc

    res_m_toA2B = rs.comp1.LIST(A_toA2B, B_toA2B, C_toA2B).Conc
            
    res_m_so2d = rs.patch1.LIST(A_so2d, B_so2d).Count

    sim.toSave(
        res_m_foi, res_m_for, res_m_soA2, res_m_soAA, res_m_soAB, 
        res_m_toA3, res_m_toA2B, res_m_so2d, 
        dt=DT
    )
    
    for i in range (0, NITER_max):
        sim.newRun()
            
        if i < NITER_foi: 
            sim.comp1.A_foi.Count = N_foi
        else:
            res_m_foi._saveWithTpnts([INT+1])
        
        if i < NITER_for: 
            sim.comp1.A_for.Count = COUNT_for
            sim.comp1.B_for.Count = 0.0
        else:
            res_m_for._saveWithTpnts([INT+1])

        if i < NITER_soA2:
            sim.comp1.A_soA2.Conc = CONCA_soA2
        else:
            res_m_soA2._saveWithTpnts([INT+1])

        if i < NITER_soAA:
            sim.comp1.A_soAA.Conc = CONCA_soAA
            sim.comp1.B_soAA.Conc = CONCB_soAA
        else:
            res_m_soAA._saveWithTpnts([INT+1])
        
        if i < NITER_soAB:
            sim.comp1.A_soAB.Conc = CONCA_soAB
            sim.comp1.B_soAB.Conc = CONCB_soAB
        else:
            res_m_soAB._saveWithTpnts([INT+1])
        
        if i < NITER_toA3:
            sim.comp1.A_toA3.Conc = CONCA_toA3
        else:
            res_m_toA3._saveWithTpnts([INT+1])

        if i < NITER_toA2B:
            sim.comp1.A_toA2B.Conc = CONCA_toA2B
            sim.comp1.B_toA2B.Conc = CONCB_toA2B
        else:
            res_m_toA2B._saveWithTpnts([INT+1])
                                    
        if i < NITER_so2d:
            sim.patch1.A_so2d.Count = COUNTA_so2d
            sim.patch1.B_so2d.Count = COUNTB_so2d
        else:
            res_m_so2d._saveWithTpnts([INT+1])            

        sim.run(INT)


    global mean_res_foi, std_res_foi
    mean_res_foi = numpy.mean(res_m_foi.data[:NITER_foi], 0)
    std_res_foi = numpy.std(res_m_foi.data, 0)

    global mean_res_for
    mean_res_for = numpy.mean(res_m_for.data[:NITER_for], 0)

    global mean_res_soA2
    mean_res_soA2 = numpy.mean(res_m_soA2.data[:NITER_soA2], 0)
    
    global mean_res_soAA
    mean_res_soAA = numpy.mean(res_m_soAA.data[:NITER_soAA], 0)

    global mean_res_soAB
    mean_res_soAB = numpy.mean(res_m_soAB.data[:NITER_soAB], 0)

    global mean_res_toA3
    mean_res_toA3 = numpy.mean(res_m_toA3.data[:NITER_toA3], 0)

    global mean_res_toA2B
    mean_res_toA2B = numpy.mean(res_m_toA2B.data[:NITER_toA2B], 0)
    
    global mean_res_so2d
    mean_res_so2d = numpy.mean(res_m_so2d.data[:NITER_so2d], 0)

    global tpnts, ntpnts
    tpnts = res_m_foi.time[0]
    ntpnts = len(tpnts)

    global ran_sim
    ran_sim = True


# Tests follow:

##################### First order irreversible #########################

def test_foi():
    "Reaction - First order, irreversible (TetOpSplit)"

    if not ran_sim:
        run_sim()

    m_tol = 0
    for i in range(ntpnts):
        if i == 0:
            continue
        analy = N_foi*math.exp(-KCST_foi*tpnts[i])
        std = math.pow((N_foi*(math.exp(-KCST_foi*tpnts[i]))*(1-(math.exp(-KCST_foi*tpnts[i])))), 0.5)
        
        assert tol_funcs.tolerable(analy, mean_res_foi[i], tolerance_foi)
        assert tol_funcs.tolerable(std, std_res_foi[i], tolerance_foi)

####################### First order reversible #########################

def test_for():
    "Reaction - First order, reversible (TetOpSplit)"

    if not ran_sim:
        run_sim()

    Aeq = COUNT_for*(KCST_b_for/KCST_f_for)/(1+(KCST_b_for/KCST_f_for))/(VOL*Avogadro * 1e3)*1e6
    Beq = (COUNT_for/(VOL*Avogadro * 1e3))*1e6 -Aeq
    for i in range(ntpnts):
        if i < 7:
            continue
        assert tol_funcs.tolerable(mean_res_for[i,0], Aeq, tolerance_for)
        assert tol_funcs.tolerable(mean_res_for[i,1], Beq, tolerance_for)

####################### Second order irreversible A2 ###################

def test_soA2():
    "Reaction - Second order, irreversible, 2A->C (TetOpSplit)"
    
    if not ran_sim:
        run_sim()
    
    invA = numpy.zeros(ntpnts)
    lineA  = numpy.zeros(ntpnts)
    for i in range(ntpnts):
        invA[i] = (1.0/mean_res_soA2[i][0])
        lineA[i] = (1.0/CONCA_soA2 +((tpnts[i]*2*KCST_soA2)))
        assert tol_funcs.tolerable(invA[i], lineA[i], tolerance_soA2)

####################### Third order irreversible A3 ###################

def test_toA3():
    "Reaction - Third order, irreversible, 3A->C (TetOpSplit)"
    
    if not ran_sim:
        run_sim()
    
    inv2A = numpy.zeros(ntpnts)
    lineA  = numpy.zeros(ntpnts)
    
    for i in range(ntpnts):
        inv2A[i] = (1.0/(mean_res_toA3[i][0]**2))
        lineA[i] = (1.0/(CONCA_toA3**2) +((tpnts[i]*6*KCST_toA3)))        
        assert tol_funcs.tolerable(inv2A[i], lineA[i], tolerance_toA3)

####################### Third order irreversible A3 ###################

def test_toA2B():
    "Reaction - Third order, irreversible, 2A+B->C (TetOpSplit)"
    
    if not ran_sim:
        run_sim()

    A0 = CONCA_toA2B
    B0 = CONCB_toA2B
    
    delta_AB = A0-2*B0
    delta_BA = 2*B0-A0
        
    kt = numpy.zeros(ntpnts)
    lineA  = numpy.zeros(ntpnts)
    for i in range(1, ntpnts):  
        A = mean_res_toA2B[i][0]
        B = mean_res_toA2B[i][1]
        lineA[i] = (-1.0/delta_AB)*( (-1.0/delta_BA)*math.log( (B/A)/(B0/A0) ) +1.0/A - 1.0/A0 )        
        kt[i] = (tpnts[i]*KCST_toA2B)
        assert tol_funcs.tolerable(kt[i], lineA[i], tolerance_toA2B)
    
####################### Second order irreversible AA ###################

def test_soAA():
    "Reaction - Second order, irreversible, A0=B0 (TetOpSplit)"

    if not ran_sim:
        run_sim()

    invA = numpy.zeros(ntpnts)
    invB = numpy.zeros(ntpnts)
    lineA  = numpy.zeros(ntpnts)
    lineB = numpy.zeros(ntpnts)
    for i in range(ntpnts):
        invA[i] = (1.0/mean_res_soAA[i][0])
        invB[i] = (1.0/mean_res_soAA[i][1])
        lineA[i] = (1.0/CONCA_soAA +((tpnts[i]*KCST_soAA)))
        lineB[i] = (1.0/CONCB_soAA + ((tpnts[i]*KCST_soAA)))
        
        assert tol_funcs.tolerable(invA[i], lineA[i], tolerance_soAA)
        assert tol_funcs.tolerable(invB[i], lineB[i], tolerance_soAA)

####################### Second order irreversible AB ###################

def test_soAB():
    "Reaction - Second order, irreversible, A0!=B0 (TetOpSplit)"

    if not ran_sim:
        run_sim()

    lnBA_soAB = numpy.zeros(ntpnts)
    lineAB_soAB = numpy.zeros(ntpnts)
    C_soAB = CONCA_soAB-CONCB_soAB
    for i in range(ntpnts):
        A_soAB = mean_res_soAB[i][0]
        B_soAB = mean_res_soAB[i][1]
        lnBA_soAB[i] = math.log(B_soAB/A_soAB)
        lineAB_soAB[i] = math.log(CONCB_soAB/CONCA_soAB) -C_soAB*KCST_soAB*tpnts[i]
        
        assert tol_funcs.tolerable(lnBA_soAB[i], lineAB_soAB[i], tolerance_soAB)

########################################################################

def test_so2d():
    "Reaction - Second-order, irreversible, 2D (TetOpSplit)"

    if not ran_sim:
        run_sim()

    lnBA_so2d = numpy.zeros(ntpnts)
    lineAB_so2d = numpy.zeros(ntpnts)

    C_so2d = COUNTA_so2d-COUNTB_so2d

    for i in range(ntpnts):
        A_so2d = mean_res_so2d[i][0]
        B_so2d = mean_res_so2d[i][1]
        lnBA_so2d[i] = math.log(B_so2d/A_so2d)
        lineAB_so2d[i] = math.log(COUNTB_so2d/COUNTA_so2d) -C_so2d*CCST_so2d*tpnts[i]
        assert tol_funcs.tolerable(lnBA_so2d[i], lineAB_so2d[i], tolerance_so2d)

########################################################################

