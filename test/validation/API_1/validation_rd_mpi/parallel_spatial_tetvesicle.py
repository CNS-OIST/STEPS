########################################################################

# Well-mixed chemical reactions.

# AIMS: to verify STEPS well-mixed stochastic solver 'TetVesicle' computes
# reaction rates correctly when many different chemical species and
# reactions are present in the model.
# Note, while these are well-mixed validations since this is a spatial
# solver any 2nd-order reacting molecles diffuse.

# This model combines all well-mixed validations, as described at:
# http://www.biomedcentral.com/content/supplementary/1752-0509-6-36-s4.pdf
# and in the equivalent individual model scripts.

########################################################################

import unittest

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.utilities.meshio as smeshio
import steps.mpi
import steps.mpi.solver as ssolv
import steps.utilities.geom_decompose as gd

import numpy
import math
import time

from . import tol_funcs

########################################################################

class TestRDMPISpatialTetvesicle(unittest.TestCase):

    def test_spatial_tetVesicle_n4(self):

        # Set up and run the simulations once, before the tests
        # analyze the results.

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

        NITER_max = 1000

        ########################################################################

        mdl  = smod.Model()
        volsys = smod.Volsys('vsys',mdl)
        surfsys = smod.Surfsys('ssys',mdl)

        # First order reversible
        A_for = smod.Spec('A_for', mdl)
        B_for = smod.Spec('B_for', mdl)
        A_for_diff = smod.Diff('A_for_diff', volsys,A_for, 0.01e-12)
        B_for_diff = smod.Diff('B_for_diff', volsys,B_for, 0.01e-12)
        R1_for = smod.Reac('R1_for', volsys, lhs = [A_for], rhs = [B_for], kcst = KCST_f_for)
        R2_for = smod.Reac('R2_for', volsys, lhs = [B_for], rhs = [A_for], kcst = KCST_b_for)

        # Second order irreversible A2
        A_soA2 = smod.Spec('A_soA2', mdl)
        C_soA2 = smod.Spec('C_soA2', mdl)
        A_soA2_diff = smod.Diff('A_soA2_diff', volsys,A_soA2, 1e-12)
        R1_soA2 = smod.Reac('R1_soA2', volsys, lhs = [A_soA2, A_soA2], rhs = [C_soA2], kcst = KCST_soA2)

        # Second order irreversible AA
        A_soAA = smod.Spec('A_soAA', mdl)
        B_soAA = smod.Spec('B_soAA', mdl)
        C_soAA = smod.Spec('C_soAA', mdl)
        A_soAA_diff = smod.Diff('A_soAA_diff', volsys,A_soAA, 0.2e-12)
        B_soAA_diff = smod.Diff('B_soAA_diff', volsys,B_soAA, 0.2e-12)
        R1_soAA = smod.Reac('R1_soAA', volsys, lhs = [A_soAA, B_soAA], rhs = [C_soAA], kcst = KCST_soAA)

        # Second order irreversible AB
        A_soAB = smod.Spec('A_soAB', mdl)
        B_soAB = smod.Spec('B_soAB', mdl)
        C_soAB = smod.Spec('C_soAB', mdl)
        A_soAB_diff = smod.Diff('A_soAB_diff', volsys,A_soAB, 0.1e-12)
        B_soAB_diff = smod.Diff('B_soAB_diff', volsys,B_soAB, 0.1e-12)
        R1_soAB = smod.Reac('R1_soAB', volsys, lhs = [A_soAB, B_soAB], rhs = [C_soAB], kcst = KCST_soAB)

        # Third order irreversible A3
        A_toA3 = smod.Spec('A_toA3', mdl)
        C_toA3 = smod.Spec('C_toA3', mdl)
        A_soA3_diff = smod.Diff('A_soA3_diff', volsys,A_toA3, 0.2e-12)
        R1_toA3 = smod.Reac('R1_toA3', volsys, lhs = [A_toA3, A_toA3, A_toA3], rhs = [C_toA3], kcst = KCST_toA3)

        # Third order irreversible A2B
        A_toA2B = smod.Spec('A_toA2B', mdl)
        B_toA2B = smod.Spec('B_toA2B', mdl)
        C_toA2B = smod.Spec('C_toA2B', mdl)
        A_soA2B_diff = smod.Diff('A_soA2B_diff', volsys,A_toA2B, 0.1e-12)
        B_soA2B_diff = smod.Diff('B_soA2B_diff', volsys,B_toA2B, 0.1e-12)
        R1_toA3 = smod.Reac('R1_toA2B', volsys, lhs = [A_toA2B, A_toA2B, B_toA2B], rhs = [C_toA2B], kcst = KCST_toA2B)

        # Second order irreversible 2D
        A_so2d = smod.Spec('A_so2d', mdl)
        B_so2d = smod.Spec('B_so2d', mdl)
        C_so2d = smod.Spec('C_so2d', mdl)
        A_so2d_diff = smod.Diff('A_so2d_diff', surfsys,A_so2d, 1.0e-12)
        B_so2d_diff = smod.Diff('B_so2d_diff', surfsys,B_so2d, 1.0e-12)
        SR1_so2d = smod.SReac('SR1_so2d', surfsys, slhs = [A_so2d, B_so2d], srhs = [C_so2d], kcst = KCST_so2d)


        mesh = smeshio.importAbaqus('validation_rd_mpi/meshes/sphere_rad1_37tets.inp', 1e-6)[0]
        VOL = mesh.getMeshVolume()

        comp1 = sgeom.TmComp('comp1', mesh, range(mesh.ntets))
        comp1.addVolsys('vsys')
        patch_tris = mesh.getSurfTris()
        patch1 = sgeom.TmPatch('patch1', mesh, patch_tris, comp1)
        patch1.addSurfsys('ssys')

        CCST_so2d = KCST_so2d/(6.02214179e23*patch1.getArea())

        rng = srng.create('r123', 512)
        rng.initialize(100)


        tet_hosts = gd.binTetsByAxis(mesh, steps.mpi.nhosts)
        tri_hosts = gd.partitionTris(mesh, tet_hosts, patch_tris)

        sim = ssolv.TetVesicle(mdl, mesh, rng, ssolv.EF_NONE)
        sim.setOutputSync(True, 0);
        sim.reset()


        global tpnts, ntpnts
        tpnts = numpy.arange(0.0, INT, DT)
        ntpnts = tpnts.shape[0]


        res_m_for = numpy.zeros([NITER_for, ntpnts, 2])

        res_m_soA2 = numpy.zeros([NITER_soA2, ntpnts, 2])

        res_m_soAA = numpy.zeros([NITER_soAA, ntpnts, 3])

        res_m_soAB = numpy.zeros([NITER_soAB, ntpnts, 3])

        res_m_toA3 = numpy.zeros([NITER_toA3, ntpnts, 2])

        res_m_toA2B = numpy.zeros([NITER_toA2B, ntpnts, 3])

        res_m_so2d = numpy.zeros([NITER_so2d, ntpnts, 3])



        for i in range (0, NITER_max):
            sim.reset()

            if i < NITER_for:
                sim.setCompSpecCount('comp1', 'A_for', COUNT_for)
                sim.setCompSpecCount('comp1', 'B_for', 0.0)

            if i < NITER_soA2:
                sim.setCompSpecConc('comp1', 'A_soA2', CONCA_soA2)

            if i < NITER_soAA:
                sim.setCompSpecConc('comp1', 'A_soAA', CONCA_soAA)
                sim.setCompSpecConc('comp1', 'B_soAA', CONCB_soAA)

            if i < NITER_soAB:
                sim.setCompSpecConc('comp1', 'A_soAB', CONCA_soAB)
                sim.setCompSpecConc('comp1', 'B_soAB', CONCB_soAB)

            if i < NITER_toA3:
                sim.setCompSpecConc('comp1', 'A_toA3', CONCA_toA3)

            if i < NITER_toA2B:
                sim.setCompSpecConc('comp1', 'A_toA2B', CONCA_toA2B)
                sim.setCompSpecConc('comp1', 'B_toA2B', CONCB_toA2B)

            if i < NITER_so2d:
                sim.setPatchSpecCount('patch1', 'A_so2d', COUNTA_so2d)
                sim.setPatchSpecCount('patch1', 'B_so2d', COUNTB_so2d)


            for t in range(0, ntpnts):
                sim.run(tpnts[t])

                if i < NITER_for:
                    res_m_for[i, t, 0] = sim.getCompSpecConc('comp1', 'A_for')*1e6
                    res_m_for[i, t, 1] = sim.getCompSpecConc('comp1', 'B_for')*1e6

                if i < NITER_soA2:
                    res_m_soA2[i, t, 0] = sim.getCompSpecConc('comp1', 'A_soA2')

                if i < NITER_soAA:
                    res_m_soAA[i, t, 0] = sim.getCompSpecConc('comp1', 'A_soAA')
                    res_m_soAA[i, t, 1] = sim.getCompSpecConc('comp1', 'B_soAA')

                if i < NITER_soAB:
                    res_m_soAB[i, t, 0] = sim.getCompSpecConc('comp1', 'A_soAB')
                    res_m_soAB[i, t, 1] = sim.getCompSpecConc('comp1', 'B_soAB')

                if i < NITER_toA3:
                    res_m_toA3[i, t, 0] = sim.getCompSpecConc('comp1', 'A_toA3')

                if i < NITER_toA2B:
                    res_m_toA2B[i, t, 0] = sim.getCompSpecConc('comp1', 'A_toA2B')
                    res_m_toA2B[i, t, 1] = sim.getCompSpecConc('comp1', 'B_toA2B')
                    res_m_toA2B[i, t, 2] = sim.getCompSpecConc('comp1', 'C_toA2B')

                if i < NITER_so2d:
                    res_m_so2d[i, t, 0] = sim.getPatchSpecCount('patch1', 'A_so2d')
                    res_m_so2d[i, t, 1] = sim.getPatchSpecCount('patch1', 'B_so2d')


        global mean_res_for
        mean_res_for = numpy.mean(res_m_for, 0)

        global mean_res_soA2
        mean_res_soA2 = numpy.mean(res_m_soA2, 0)

        global mean_res_soAA
        mean_res_soAA = numpy.mean(res_m_soAA, 0)

        global mean_res_soAB
        mean_res_soAB = numpy.mean(res_m_soAB, 0)

        global mean_res_toA3
        mean_res_toA3 = numpy.mean(res_m_toA3, 0)

        global mean_res_toA2B
        mean_res_toA2B = numpy.mean(res_m_toA2B, 0)

        global mean_res_so2d
        mean_res_so2d = numpy.mean(res_m_so2d, 0)

        global ran_sim
        ran_sim = True

    # Tests follow:

    ####################### First order reversible #########################

        "Reaction - First order, reversible (TetVesicle)"

        Aeq = COUNT_for*(KCST_b_for/KCST_f_for)/(1+(KCST_b_for/KCST_f_for))/(VOL*6.0221415e26)*1e6
        Beq = (COUNT_for/(VOL*6.0221415e26))*1e6 -Aeq
        for i in range(ntpnts):
            if i < 7: continue
            self.assertTrue(tol_funcs.tolerable(mean_res_for[i,0], Aeq, tolerance_for))
            self.assertTrue(tol_funcs.tolerable(mean_res_for[i,1], Beq, tolerance_for))

    ####################### Second order irreversible A2 ###################

        "Reaction - Second order, irreversible, 2A->C (TetVesicle)"

        invA = numpy.zeros(ntpnts)
        lineA  = numpy.zeros(ntpnts)
        for i in range(ntpnts):
            invA[i] = (1.0/mean_res_soA2[i][0])
            lineA[i] = (1.0/CONCA_soA2 +((tpnts[i]*2*KCST_soA2)))
            self.assertTrue(tol_funcs.tolerable(invA[i], lineA[i], tolerance_soA2))

    ####################### Third order irreversible A3 ###################

        "Reaction - Third order, irreversible, 3A->C (TetVesicle)"

        inv2A = numpy.zeros(ntpnts)
        lineA  = numpy.zeros(ntpnts)

        for i in range(ntpnts):
            inv2A[i] = (1.0/(mean_res_toA3[i][0]**2))
            lineA[i] = (1.0/(CONCA_toA3**2) +((tpnts[i]*6*KCST_toA3)))
            self.assertTrue(tol_funcs.tolerable(inv2A[i], lineA[i], tolerance_toA3))

    ####################### Third order irreversible A3 ###################

        "Reaction - Third order, irreversible, 2A+B->C (TetVesicle)"

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
            self.assertTrue(tol_funcs.tolerable(kt[i], lineA[i], tolerance_toA2B))

    ####################### Second order irreversible AA ###################

        "Reaction - Second order, irreversible, A0=B0 (TetVesicle)"

        invA = numpy.zeros(ntpnts)
        invB = numpy.zeros(ntpnts)
        lineA  = numpy.zeros(ntpnts)
        lineB = numpy.zeros(ntpnts)
        for i in range(ntpnts):
            invA[i] = (1.0/mean_res_soAA[i][0])
            invB[i] = (1.0/mean_res_soAA[i][1])
            lineA[i] = (1.0/CONCA_soAA +((tpnts[i]*KCST_soAA)))
            lineB[i] = (1.0/CONCB_soAA + ((tpnts[i]*KCST_soAA)))

            self.assertTrue(tol_funcs.tolerable(invA[i], lineA[i], tolerance_soAA))
            self.assertTrue(tol_funcs.tolerable(invB[i], lineB[i], tolerance_soAA))

    ####################### Second order irreversible AB ###################

        "Reaction - Second order, irreversible, A0!=B0 (TetVesicle)"

        lnBA_soAB = numpy.zeros(ntpnts)
        lineAB_soAB = numpy.zeros(ntpnts)
        C_soAB = CONCA_soAB-CONCB_soAB
        for i in range(ntpnts):
            A_soAB = mean_res_soAB[i][0]
            B_soAB = mean_res_soAB[i][1]
            lnBA_soAB[i] = math.log(B_soAB/A_soAB)
            lineAB_soAB[i] = math.log(CONCB_soAB/CONCA_soAB) -C_soAB*KCST_soAB*tpnts[i]

            self.assertTrue(tol_funcs.tolerable(lnBA_soAB[i], lineAB_soAB[i], tolerance_soAB))

    ########################################################################

        "Reaction - Second-order, irreversible, 2D (TetVesicle)"

        lnBA_so2d = numpy.zeros(ntpnts)
        lineAB_so2d = numpy.zeros(ntpnts)

        C_so2d = COUNTA_so2d-COUNTB_so2d

        for i in range(ntpnts):
            A_so2d = mean_res_so2d[i][0]
            B_so2d = mean_res_so2d[i][1]
            lnBA_so2d[i] = math.log(B_so2d/A_so2d)
            lineAB_so2d[i] = math.log(COUNTB_so2d/COUNTA_so2d) -C_so2d*CCST_so2d*tpnts[i]
            self.assertTrue(tol_funcs.tolerable(lnBA_so2d[i], lineAB_so2d[i], tolerance_so2d))

    ########################################################################

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestRDMPISpatialTetvesicle))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
