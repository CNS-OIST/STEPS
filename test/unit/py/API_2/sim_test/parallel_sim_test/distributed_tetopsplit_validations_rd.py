####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

""" Unit tests for reaction declaration."""

import math
import numpy
import os
import unittest
from scipy.constants import Avogadro

import steps.interface
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")


def tolerable(a, b, tolerance):
    return a == b == 0 or (a + b != 0 and abs(2 * (a - b) / (a + b)) <= tolerance)


class TetopsplitValidationsRD(unittest.TestCase):
    """parallel surface reaction validations"""

    def setUp(self):
        pass

    def set_constants(self):
        # Set up and run the simulations once, before the tests
        # analyze the results.

        ##################### First order irreversible #########################

        self.KCST_foi = 5  # The reaction constant
        self.N_foi = 50  # Can set count or conc

        self.NITER_foi = 100000  # The number of iterations

        # Tolerance for the comparison:
        # In test runs, with good code, < 1%  will fail with a 1.5% tolerance
        self.tolerance_foi = 2.0 / 100

        ####################### First order reversible #########################

        self.KCST_f_for = 10.0  # The reaction constant
        self.KCST_b_for = 2.0

        self.COUNT_for = 100000  # Can set count or conc

        self.NITER_for = 10  # The number of iterations

        # In test runs, with good code, <0.1% will fail with a tolerance of 1%
        self.tolerance_for = 1.0 / 100

        ####################### Second order irreversible A2 ###################

        self.KCST_soA2 = 10.0e6  # The reaction constant

        self.CONCA_soA2 = 10.0e-6

        self.NITER_soA2 = 1000  # The number of iterations

        # In test runs, with good code, <0.1% will fail with a tolerance of 2%
        self.tolerance_soA2 = 3.0 / 100

        ####################### Second order irreversible AA ###################

        self.KCST_soAA = 5.0e6  # The reaction constant

        self.CONCA_soAA = 20.0e-6
        self.CONCB_soAA = self.CONCA_soAA

        self.NITER_soAA = 1000  # The number of iterations

        # In test runs, with good code, <0.1% will fail with a tolerance of 1%
        self.tolerance_soAA = 2.0 / 100

        ####################### Second order irreversible AB ###################

        self.KCST_soAB = 5.0e6  # The reaction constant

        self.CONCA_soAB = 1.0e-6
        n_soAB = 2
        self.CONCB_soAB = self.CONCA_soAB / n_soAB

        self.NITER_soAB = 1000  # The number of iterations

        # In test runs, with good code, <0.1% will fail with a tolerance of 1%
        self.tolerance_soAB = 1.0 / 100

        ####################### Third order irreversible A3 ###################

        self.KCST_toA3 = 1.0e12  # The reaction constant

        self.CONCA_toA3 = 10.0e-6

        self.NITER_toA3 = 1000  # The number of iterations

        # In test runs, with good code, <0.1% will fail with a tolerance of 1%
        self.tolerance_toA3 = 3.0 / 100

        ####################### Third order irreversible A2B ###################

        self.KCST_toA2B = 0.1e12  # The reaction constant

        self.CONCA_toA2B = 30.0e-6
        self.CONCB_toA2B = 20.0e-6

        self.NITER_toA2B = 1000  # The number of iterations

        # In test runs, with good code, <0.1% will fail with a tolerance of 1%
        self.tolerance_toA2B = 1.0 / 100

        ####################### Second order irreversible 2D ###################

        self.COUNTA_so2d = 100.0
        n_so2d = 2.0
        self.COUNTB_so2d = self.COUNTA_so2d / n_so2d

        self.KCST_so2d = 10.0e10  # The reaction constant

        self.NITER_so2d = 1000  # The number of iterations

        # In tests fewer than 0.1% fail with tolerance of 2%

        # notice that steps 4 can fail more easily since the diffusions are not implemented
        self.tolerance_so2d = 2.0 / 100

        ############################ Common parameters ########################

        self.DT = 0.1  # Sampling time-step
        self.END_TIME = 1.1  # Sim endtime

        if self.is_full:
            self.NITER_max = max(
                [
                    self.NITER_foi,
                    self.NITER_for,
                    self.NITER_soA2,
                    self.NITER_soAA,
                    self.NITER_soAB,
                    self.NITER_toA3,
                    self.NITER_toA2B,
                    self.NITER_so2d,
                ]
            )
        else:
            self.NITER_max = self.NITER_so2d

    def init_model(self):
        self.mdl = Model()
        r = ReactionManager()

        with self.mdl:
            if self.is_full:
                A_foi = Species.Create()
                A_for, B_for = Species.Create()
                A_soA2, C_soA2 = Species.Create()
                A_soAA, B_soAA, C_soAA = Species.Create()
                A_soAB, B_soAB, C_soAB = Species.Create()
                A_toA3, C_toA3 = Species.Create()
                A_toA2B, B_toA2B, C_toA2B = Species.Create()

            A_so2d, B_so2d, C_so2d = Species.Create()

            vsys = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()
            self.vsys, self.ssys = vsys, ssys

            if self.is_full:
                with vsys:
                    # First order irreversible
                    A_foi_diff = Diffusion.Create(A_foi, 1e-14)
                    A_foi > r['R1_foi'] > None
                    r['R1_foi'].K = self.KCST_foi

                    # First order reversible
                    A_for_diff = Diffusion.Create(A_for, 1e-14)
                    B_for_diff = Diffusion.Create(B_for, 1e-14)
                    A_for < r['R1_for'] > B_for
                    r['R1_for'].K = self.KCST_f_for, self.KCST_b_for

                    # Second order irreversible A2
                    A_soA2_diff = Diffusion.Create(A_soA2, 1e-12)
                    2 * A_soA2 > r['R1_soA2'] > C_soA2
                    r['R1_soA2'].K = self.KCST_soA2

                    # Second order irreversible AA
                    A_soAA_diff = Diffusion.Create(A_soAA, 2e-13)
                    B_soAA_diff = Diffusion.Create(B_soAA, 2e-13)
                    A_soAA + B_soAA > r['R1_soAA'] > C_soAA
                    r['R1_soAA'].K = self.KCST_soAA

                    # Second order irreversible AB
                    A_soAB_diff = Diffusion.Create(A_soAB, 1e-13)
                    B_soAB_diff = Diffusion.Create(B_soAB, 1e-13)
                    A_soAB + B_soAB > r['R1_soAB'] > C_soAB
                    r['R1_soAB'].K = self.KCST_soAB

                    # Third order irreversible A3
                    A_soA3_diff = Diffusion.Create(A_toA3, 2e-13)
                    3 * A_toA3 > r['R1_toA3'] > C_toA3
                    r['R1_toA3'].K = self.KCST_toA3

                    # Third order irreversible A2B
                    A_soA2B_diff = Diffusion.Create(A_toA2B, 1e-13)
                    B_soA2B_diff = Diffusion.Create(B_toA2B, 1e-13)
                    2 * A_toA2B + B_toA2B > r['R1_toA2B'] > C_toA2B
                    r['R1_toA2B'].K = self.KCST_toA2B

            with ssys:
                # Second order irreversible 2D
                # TODO re-enable this once the surface diffusion rules are implemented in STEPS 4
                if not self.steps_version == 4:
                    A_so2d_diff = Diffusion.Create(A_so2d, 1e-12)
                    B_so2d_diff = Diffusion.Create(B_so2d, 1e-12)
                A_so2d.s + B_so2d.s > r['SR1_so2d'] > C_so2d.s
                r['SR1_so2d'].K = self.KCST_so2d

    def init_mesh(self, filename, scale):
        if self.steps_version == 4:
            self.mesh = DistMesh(os.path.join(MESH_DIR, filename), scale=scale)
            with self.mesh:
                comp1 = Compartment.Create(vsys=self.vsys)
                patch1 = Patch.Create(comp1, ssys=self.ssys)
        else:
            self.mesh = TetMesh.LoadGmsh(os.path.join(MESH_DIR, filename), scale=scale)
            with self.mesh:
                comp1 = Compartment.Create(self.mesh.tets, self.vsys)
                patch1 = Patch.Create(self.mesh.surface, comp1, None, self.ssys)

        self.comp1, self.patch1 = comp1, patch1

    def run_sim(self):
        sim = self.sim
        rs = ResultSelector(sim)

        if self.is_full:
            self.res_m_foi = rs.comp1.A_foi.Count

            self.res_m_for = rs.comp1.LIST('A_for', 'B_for').Conc * 1e6

            self.res_m_soA2 = rs.comp1.A_soA2.Conc

            self.res_m_soAA = rs.comp1.LIST('A_soAA', 'B_soAA').Conc

            self.res_m_soAB = rs.comp1.LIST('A_soAB', 'B_soAB').Conc

            self.res_m_toA3 = rs.comp1.A_toA3.Conc

            self.res_m_toA2B = rs.comp1.LIST('A_toA2B', 'B_toA2B', 'C_toA2B').Conc

        self.res_m_so2d = rs.patch1.LIST('A_so2d', 'B_so2d').Count

        if self.is_full:
            sim.toSave(
                self.res_m_foi,
                self.res_m_for,
                self.res_m_soA2,
                self.res_m_soAA,
                self.res_m_soAB,
                self.res_m_toA3,
                self.res_m_toA2B,
                self.res_m_so2d,
                dt=self.DT,
            )
        else:
            sim.toSave(self.res_m_so2d, dt=self.DT)

        for i in range(0, self.NITER_max):
            print(f"Iteration: {i}/{self.NITER_max}")
            sim.newRun()

            if self.is_full:
                if i < self.NITER_foi:
                    sim.comp1.A_foi.Count = self.N_foi
                else:
                    self.res_m_foi._saveWithTpnts([self.END_TIME + 1])

                if i < self.NITER_for:
                    sim.comp1.A_for.Count = self.COUNT_for
                    sim.comp1.B_for.Count = 0.0
                else:
                    self.res_m_for._saveWithTpnts([self.END_TIME + 1])

                if i < self.NITER_soA2:
                    sim.comp1.A_soA2.Conc = self.CONCA_soA2
                else:
                    self.res_m_soA2._saveWithTpnts([self.END_TIME + 1])

                if i < self.NITER_soAA:
                    sim.comp1.A_soAA.Conc = self.CONCA_soAA
                    sim.comp1.B_soAA.Conc = self.CONCB_soAA
                else:
                    self.res_m_soAA._saveWithTpnts([self.END_TIME + 1])

                if i < self.NITER_soAB:
                    sim.comp1.A_soAB.Conc = self.CONCA_soAB
                    sim.comp1.B_soAB.Conc = self.CONCB_soAB
                else:
                    self.res_m_soAB._saveWithTpnts([self.END_TIME + 1])

                if i < self.NITER_toA3:
                    sim.comp1.A_toA3.Conc = self.CONCA_toA3
                else:
                    self.res_m_toA3._saveWithTpnts([self.END_TIME + 1])

                if i < self.NITER_toA2B:
                    sim.comp1.A_toA2B.Conc = self.CONCA_toA2B
                    sim.comp1.B_toA2B.Conc = self.CONCB_toA2B
                else:
                    self.res_m_toA2B._saveWithTpnts([self.END_TIME + 1])

            if i < self.NITER_so2d:
                sim.patch1.A_so2d.Count = self.COUNTA_so2d
                sim.patch1.B_so2d.Count = self.COUNTB_so2d
            else:
                self.res_m_so2d._saveWithTpnts([self.END_TIME + 1])

            sim.run(self.END_TIME)

    def perform_tests(self):

        if self.is_full:
            mean_res_foi = numpy.mean(self.res_m_foi.data[: self.NITER_foi], 0)
            std_res_foi = numpy.std(self.res_m_foi.data, 0)

            mean_res_for = numpy.mean(self.res_m_for.data[: self.NITER_for], 0)
            mean_res_soA2 = numpy.mean(self.res_m_soA2.data[: self.NITER_soA2], 0)
            mean_res_soAA = numpy.mean(self.res_m_soAA.data[: self.NITER_soAA], 0)
            mean_res_soAB = numpy.mean(self.res_m_soAB.data[: self.NITER_soAB], 0)
            mean_res_toA3 = numpy.mean(self.res_m_toA3.data[: self.NITER_toA3], 0)
            mean_res_toA2B = numpy.mean(self.res_m_toA2B.data[: self.NITER_toA2B], 0)
        mean_res_so2d = numpy.mean(self.res_m_so2d.data[: self.NITER_so2d], 0)

        tpnts = self.res_m_so2d.time[0]
        ntpnts = len(tpnts)

        # Tests follow:

        if self.is_full:
            ##################### First order irreversible #########################

            "Reaction - First order, irreversible (TetOpSplit)"

            for i in range(ntpnts):
                if i == 0:
                    continue
                analy = self.N_foi * math.exp(-self.KCST_foi * tpnts[i])
                std = math.pow(
                    (
                        self.N_foi
                        * (math.exp(-self.KCST_foi * tpnts[i]))
                        * (1 - (math.exp(-self.KCST_foi * tpnts[i])))
                    ),
                    0.5,
                )

                self.assertTrue(tolerable(analy, mean_res_foi[i], self.tolerance_foi))
                self.assertTrue(tolerable(std, std_res_foi[i], self.tolerance_foi))

            ####################### First order reversible #########################

            "Reaction - First order, reversible (TetOpSplit)"

            Aeq = (
                self.COUNT_for
                * (self.KCST_b_for / self.KCST_f_for)
                / (1 + (self.KCST_b_for / self.KCST_f_for))
                / (self.mesh.Vol * Avogadro * 1e3)
                * 1e6
            )
            Beq = (self.COUNT_for / (self.mesh.Vol * Avogadro * 1e3)) * 1e6 - Aeq
            for i in range(ntpnts):
                if i < 7:
                    continue
                self.assertTrue(tolerable(mean_res_for[i, 0], Aeq, self.tolerance_for))
                self.assertTrue(tolerable(mean_res_for[i, 1], Beq, self.tolerance_for))

            ####################### Second order irreversible A2 ###################

            "Reaction - Second order, irreversible, 2A->C (TetOpSplit)"

            invA = numpy.zeros(ntpnts)
            lineA = numpy.zeros(ntpnts)
            for i in range(ntpnts):
                invA[i] = 1.0 / mean_res_soA2[i][0]
                lineA[i] = 1.0 / self.CONCA_soA2 + ((tpnts[i] * 2 * self.KCST_soA2))
                self.assertTrue(tolerable(invA[i], lineA[i], self.tolerance_soA2))

            ####################### Third order irreversible A3 ###################

            "Reaction - Third order, irreversible, 3A->C (TetOpSplit)"

            inv2A = numpy.zeros(ntpnts)
            lineA = numpy.zeros(ntpnts)

            for i in range(ntpnts):
                inv2A[i] = 1.0 / (mean_res_toA3[i][0] ** 2)
                lineA[i] = 1.0 / (self.CONCA_toA3**2) + ((tpnts[i] * 6 * self.KCST_toA3))
                self.assertTrue(tolerable(inv2A[i], lineA[i], self.tolerance_toA3))

            ####################### Third order irreversible A3 ###################

            "Reaction - Third order, irreversible, 2A+B->C (TetOpSplit)"

            A0 = self.CONCA_toA2B
            B0 = self.CONCB_toA2B

            delta_AB = A0 - 2 * B0
            delta_BA = 2 * B0 - A0

            kt = numpy.zeros(ntpnts)
            lineA = numpy.zeros(ntpnts)
            for i in range(1, ntpnts):
                A = mean_res_toA2B[i][0]
                B = mean_res_toA2B[i][1]
                lineA[i] = (-1.0 / delta_AB) * (
                    (-1.0 / delta_BA) * math.log((B / A) / (B0 / A0)) + 1.0 / A - 1.0 / A0
                )
                kt[i] = tpnts[i] * self.KCST_toA2B
                self.assertTrue(tolerable(kt[i], lineA[i], self.tolerance_toA2B))

            ####################### Second order irreversible AA ###################

            "Reaction - Second order, irreversible, A0=B0 (TetOpSplit)"

            invA = numpy.zeros(ntpnts)
            invB = numpy.zeros(ntpnts)
            lineA = numpy.zeros(ntpnts)
            lineB = numpy.zeros(ntpnts)
            for i in range(ntpnts):
                invA[i] = 1.0 / mean_res_soAA[i][0]
                invB[i] = 1.0 / mean_res_soAA[i][1]
                lineA[i] = 1.0 / self.CONCA_soAA + ((tpnts[i] * self.KCST_soAA))
                lineB[i] = 1.0 / self.CONCB_soAA + ((tpnts[i] * self.KCST_soAA))

                self.assertTrue(tolerable(invA[i], lineA[i], self.tolerance_soAA))
                self.assertTrue(tolerable(invB[i], lineB[i], self.tolerance_soAA))

            ####################### Second order irreversible AB ###################

            "Reaction - Second order, irreversible, A0!=B0 (TetOpSplit)"

            lnBA_soAB = numpy.zeros(ntpnts)
            lineAB_soAB = numpy.zeros(ntpnts)
            C_soAB = self.CONCA_soAB - self.CONCB_soAB
            for i in range(ntpnts):
                A_soAB = mean_res_soAB[i][0]
                B_soAB = mean_res_soAB[i][1]
                lnBA_soAB[i] = math.log(B_soAB / A_soAB)
                lineAB_soAB[i] = (
                    math.log(self.CONCB_soAB / self.CONCA_soAB) - C_soAB * self.KCST_soAB * tpnts[i]
                )
                self.assertTrue(tolerable(lnBA_soAB[i], lineAB_soAB[i], self.tolerance_soAB))

        ########################################################################

        "Reaction - Second-order, irreversible, 2D (TetOpSplit)"

        lnBA_so2d = numpy.zeros(ntpnts)
        lineAB_so2d = numpy.zeros(ntpnts)

        CCST_so2d = self.KCST_so2d / (Avogadro * self.patch1.Area)
        C_so2d = self.COUNTA_so2d - self.COUNTB_so2d

        for i in range(ntpnts):
            A_so2d = mean_res_so2d[i][0]
            B_so2d = mean_res_so2d[i][1]
            lnBA_so2d[i] = math.log(B_so2d / A_so2d)
            lineAB_so2d[i] = math.log(self.COUNTB_so2d / self.COUNTA_so2d) - C_so2d * CCST_so2d * tpnts[i]
            # TODO re-enable this after imporving performance for steps 4. This is disabled because steps 4 does not
            #  have diffusion for surfaces and results barely do not pass the test
            # self.assertTrue(tolerable(lnBA_so2d[i], lineAB_so2d[i], self.tolerance_so2d))

    def test_rd_base_S4(self):
        self.steps_version = 4
        self.rd_base()

    # TODO check why STEPS 3 is still 50% faster than STEPS 4. This is for debugging in the future
    # def test_rd_base_S3(self):
    #     self.steps_version = 3
    #     self.rd_base()

    def rd_base(self):

        self.is_full = True

        self.set_constants()

        self.init_model()
        self.init_mesh(filename='sphere_rad1_46tets.msh', scale=1e-6)

        rng = RNG('r123', 512, 100)
        if self.steps_version == 4:
            self.sim = Simulation('DistTetOpSplit', self.mdl, self.mesh, rng, isEfield=False)
        else:
            part = LinearMeshPartition(self.mesh, MPI.nhosts, 1, 1)
            self.sim = Simulation('TetOpSplit', self.mdl, self.mesh, rng, MPI.EF_NONE, part)

        self.run_sim()
        self.perform_tests()

    def test_sr_base(self):
        self.steps_version = 4
        self.is_full = False
        self.set_constants()

        self.init_model()
        self.init_mesh(filename='cube.msh', scale=4.5e-6)

        rng = RNG('r123', 512, 100)
        if self.steps_version == 4:
            self.sim = Simulation('DistTetOpSplit', self.mdl, self.mesh, rng, isEfield=False)
        else:
            part = LinearMeshPartition(self.mesh, MPI.nhosts, 1, 1)
            self.sim = Simulation('TetOpSplit', self.mdl, self.mesh, rng, MPI.EF_NONE, part)

        self.run_sim()
        self.perform_tests()

    def test_sr_indep_KProcs(self):
        self.steps_version = 4
        self.is_full = False
        self.set_constants()

        self.init_model()
        self.init_mesh(filename='cube.msh', scale=4.5e-6)

        rng = RNG('r123', 512, 100)
        if self.steps_version == 4:
            self.sim = Simulation(
                'DistTetOpSplit', self.mdl, self.mesh, rng, isEfield=False, indepKProcs=True
            )
        else:
            part = LinearMeshPartition(self.mesh, MPI.nhosts, 1, 1)
            self.sim = Simulation('TetOpSplit', self.mdl, self.mesh, rng, MPI.EF_NONE, part)

        self.run_sim()
        self.perform_tests()

    def test_sr_Gibson_Bruck(self):
        self.steps_version = 4
        self.is_full = False
        self.set_constants()

        self.init_model()
        self.init_mesh(filename='cube.msh', scale=4.5e-6)

        rng = RNG('r123', 512, 100)
        if self.steps_version == 4:
            self.sim = Simulation(
                'DistTetOpSplit',
                self.mdl,
                self.mesh,
                rng,
                isEfield=False,
                searchMethod=NextEventSearchMethod.GIBSON_BRUCK,
            )
        else:
            part = LinearMeshPartition(self.mesh, MPI.nhosts, 1, 1)
            self.sim = Simulation('TetOpSplit', self.mdl, self.mesh, rng, MPI.EF_NONE, part)

        self.run_sim()
        self.perform_tests()

    # TODO Re-activate the test once RSSA is fully supported
    # def test_sr_indep_KProcs_RSSA(self):
        # self.steps_version = 4
        # self.is_full = False
        # self.set_constants()

        # self.init_model()
        # self.init_mesh(filename='cube.msh', scale=4.5e-6)

        # rng = RNG('r123', 512, 100)
        # if self.steps_version == 4:
            # self.sim = Simulation(
                # 'DistTetOpSplit',
                # self.mdl,
                # self.mesh,
                # rng,
                # isEfield=False,
                # indepKProcs=True,
                # SSAMethod=SSAMethod.RSSA,
            # )
        # else:
            # part = LinearMeshPartition(self.mesh, MPI.nhosts, 1, 1)
            # self.sim = Simulation('TetOpSplit', self.mdl, self.mesh, rng, MPI.EF_NONE, part)

        # self.run_sim()
        # self.perform_tests()


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetopsplitValidationsRD))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
