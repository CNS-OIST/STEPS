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

"""
Integration test inspired by the caburst test.

- reactions
- diffusions
- GHK reactions
- multi compartment
- efield
- petsc keys
 """

import os
import unittest
import numpy as np
import scipy.constants as spc

import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")

N_AVOGADRO = spc.physical_constants['Avogadro constant'][0]


class DistTetopsplitMiniCaburst(unittest.TestCase):
    """Test loosely based on ca burst"""

    def setConstants(self, trash_eqv=False):
        self.DT = 1e-4
        self.EFIELD_DT = 1e-5
        self.NTSTEPS = 10
        self.END_TIME = self.DT * self.NTSTEPS

        self.petsc_options = "-ksp_rtol 1e-8 -ksp_norm_type unpreconditioned"

        self.DEFAULT_EQV = -65.0e-3
        self.TRASHED_EQV = -65.0e-2
        if trash_eqv:
            self.eqV = self.TRASHED_EQV
        else:
            self.eqV = self.DEFAULT_EQV
        self.base_conc = 1e21 / N_AVOGADRO
        self.base_temp = 273 + 30.0

    def setUpModel(self, **kwargs):
        """Set up model with reactions and diffusions"""
        self.model = Model()
        r = ReactionManager()

        with self.model:
            vsys1 = VolumeSystem.Create()
            vsys2 = VolumeSystem.Create()
            SA, SB = Species.Create()
            SC = Species.Create(valence=2)

            with vsys1:
                SC > r['R1'] > SA
                r['R1'].K = 1e3
                SA > r['R2'] > SB
                r['R2'].K = 1e3
                SB > r['R3'] > SC
                r['R3'].K = 1e3
                DC1 = Diffusion.Create(SC, 1e-7)

            with vsys2:
                SC > r['R4'] > None
                r['R4'].K = 1e3
                DC2 = Diffusion.Create(SC, 1e-7)

            ssys1 = SurfaceSystem.Create()
            Cchan1state0 = SubUnitState.Create()
            Cchan1 = Channel.Create([Cchan1state0])
            Lchan1state0 = SubUnitState.Create()
            Lchan1 = Channel.Create([Lchan1state0])
            with ssys1:
                CGHKcurr1 = GHKCurr.Create(
                    Cchan1[Cchan1state0], SC, 1.0e-20, computeflux=True, virtual_oconc=10 * self.base_conc
                )
                LOHMcurr1 = OhmicCurr.Create(Lchan1[Lchan1state0], 1e-9, self.eqV)
                SC.s > r['R5'] > SC.i
                r['R5'].K = 1e2

                SC.i > r['R6'] > SC.s

                def vRateFuncR6(V):
                    y1 = 1e1
                    y2 = 1e2
                    x1 = -0.065
                    x2 = -0.0645
                    a = (y1 - y2) / (x1 - x2)
                    q = y1 - a * x1
                    return a * V + q

                r['R6'].K = VDepRate(vRateFuncR6)

            ssysInBetween = SurfaceSystem.Create()
            CchanInBetweenState0 = SubUnitState.Create()
            CchanInBetween = Channel.Create([CchanInBetweenState0])
            with ssysInBetween:
                CGHKcurrInBetween = GHKCurr.Create(
                    CchanInBetween[CchanInBetweenState0], SC, 1.0e-20, computeflux=True
                )

        self.vsys1, self.vsys2, self.ssys1, self.ssysInBetween, self.Lchan1state0 = vsys1, vsys2, ssys1, ssysInBetween, Lchan1state0

    def setUpMeshes(self, path):
        """Create mesh and assign to it the VolumeSystem"""
        self.mesh = DistMesh(os.path.join(MESH_DIR, path), scale=1e-6)
        with self.mesh:
            comp1 = Compartment.Create(vsys=self.vsys1, conductivity=1)
            comp2 = Compartment.Create(vsys=self.vsys2, conductivity=1)
            patch1 = Patch.Create(inner=comp1, ssys=self.ssys1)
            patchInBetween = Patch.Create(inner=comp1, outer=comp2, ssys=self.ssysInBetween)
            memb1 = Membrane.Create([patch1], capacitance=1)
            membInBetween = Membrane.Create([patchInBetween], capacitance=1)

    def setUpSimulation(self):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, MPI.rank)
        self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng, isEfield=True, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)

    def setUpInitialConditions(self, trash_eqv=False):
        """Set initial conditions"""
        self.sim.Temp = self.base_temp
        self.sim.ALL(Membrane).Potential = self.eqV
        self.sim.EfieldDT = self.EFIELD_DT

        self.sim.setPetscOptions(self.petsc_options)

        self.sim.comp1.SC.Conc = self.base_conc
        self.sim.comp2.SC.Conc = self.base_conc / 2

        if trash_eqv:
            self._testModifiersOfOhmicCurrentReversalPotential()

        self.sim.patch1.Cchan1[self.model.Cchan1state0].Count = 10
        self.sim.patch1.Lchan1[self.model.Lchan1state0].Count = 1
        self.sim.patch1.SC.Count = 0
        self.sim.patchInBetween.CchanInBetween[self.model.CchanInBetweenState0].Count = 10

    def _testModifiersOfOhmicCurrentReversalPotential(self):
      """test the get/setOhmicErev methods of STEPS 4 solver"""
      tris = self.sim.TRIS(self.mesh.patch1.tris)
      ohmic_current = tris.LOHMcurr1[self.Lchan1state0]
      # test getter
      self.assertEqual(ohmic_current.Erev, self.TRASHED_EQV)
      # test modifier
      ohmic_current.Erev = self.DEFAULT_EQV
      self.assertEqual(ohmic_current.Erev, self.DEFAULT_EQV)

      self.sim.ALL(Membrane).Potential = self.DEFAULT_EQV

    def run_sim(self, **kwargs):
        """Test after ENDTIME that the counts match with the results of a previous run taken as reference"""

        self.setUpSimulation()
        tris = self.mesh.patch1.tris
        rs = ResultSelector(self.sim)
        counts = rs.comp1.LIST('SA', 'SB', 'SC').Count << rs.patch1.SC.Count
        voltages = rs.MIN(rs.VERTS(tris.verts).V) << rs.MAX(rs.VERTS(tris.verts).V)
        GHKCurr = rs.TRIS(tris).CGHKcurr1.I

        self.sim.toSave(counts, voltages, GHKCurr, dt=self.DT)
        self.sim.newRun()
        self.setUpInitialConditions(**kwargs)

        self.sim.run(self.END_TIME)

        if MPI.rank == 0:
            print("Counts: SA_comp1, SB_comp1, SC_comp1, SC_patch1\n", counts.data[0])
            print("Potential on verts: min_patch1, max_patch1\n", voltages.data[0])
            print("GHKcurr_patch1:\n", GHKCurr.data[0])
            self.assertTrue(np.isclose(voltages.data[0,-1,1], -0.0639663, rtol=0, atol=1e-4))

            # test dumpDepGraphToFile
            dep_graph_path = "dep_graph.dot"
            self.assertFalse(os.path.exists(dep_graph_path))
            self.sim.dumpDepGraphToFile("dep_graph.dot")
            self.assertTrue(os.path.exists(dep_graph_path))
            os.remove(dep_graph_path)


    def _test_tetopsplit_mini_caburst_subtest(self, mesh_file, **kwargs):
        self.setConstants(**kwargs)
        self.setUpModel(**kwargs)
        self.setUpMeshes(mesh_file)
        self.run_sim(**kwargs)

    def _test_tetopsplit_mini_caburst(self, mesh_file):
        subtests = [
            # straightforward caburst
            dict(mesh_file=mesh_file, trash_eqv=False),

            # originally set an invalid reversal potential for ohmic current
            # and fix it with get/setOhmicErev methods
            dict(mesh_file=mesh_file, trash_eqv=True),
        ]
        for subtest in subtests:
            with self.subTest(**subtest):
                self._test_tetopsplit_mini_caburst_subtest(**subtest)

    def test_tetopsplit_mini_caburst_n1(self):
        self._test_tetopsplit_mini_caburst("3tets_2patches_2comp.msh")

    def test_tetopsplit_mini_caburst_n2(self):
        self._test_tetopsplit_mini_caburst("3tets_2patches_2comp_split2/3tets_2patches_2comp")


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistTetopsplitMiniCaburst))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
