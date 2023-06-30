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

""" Unit tests for multiple compartment reaction/diffusions"""

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


def tolerable(a, b, tolerance):
    return a == b == 0 or (a + b != 0 and abs(2 * (a - b) / (a + b)) <= tolerance)


class DistTetopsplitMultiCompRd(unittest.TestCase):
    """Test multi compartments reaction-diffusions with split and whole meshes"""

    def setUp(self):
        self.setConstants()
        self.setUpModel()

    def setConstants(self):
        self.END_TIME = 1

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        self.model = Model()
        r = ReactionManager()

        with self.model:
            vsys1 = VolumeSystem.Create()
            vsys2 = VolumeSystem.Create()
            SA, SB, SC = Species.Create()
            with vsys1:
                SA + SB < r['R1'] > SC
                r['R1'].K = 1e12, 1e2
                DA1 = Diffusion.Create(SA, 1e-8)
                DB1 = Diffusion.Create(SB, 1e-8)
                DC1 = Diffusion.Create(SC, 1e-8)
            with vsys2:
                SA + SB < r['R2'] > SC
                r['R2'].K = 1e12, 1e2
                DA2 = Diffusion.Create(SA, 1e-8)
                DB2 = Diffusion.Create(SB, 1e-8)
                DC2 = Diffusion.Create(SC, 1e-8)

        self.vsys1, self.vsys2 = vsys1, vsys2

    def setUpMeshes(self, path):
        """Create mesh and assign to it the VolumeSystem"""
        # "3tets_2patches_2comp_split2/3tets_2patches_2comp"
        self.mesh = DistMesh(os.path.join(MESH_DIR, path), scale=1e-6)
        with self.mesh:
            comp1 = Compartment.Create(vsys=self.vsys1)
            comp2 = Compartment.Create(vsys=self.vsys2)

    def setUpSimulation(self):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, 7233)
        self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng, isEfield=False)

    def setUpInitialConditions(self):
        """Set initial conditions"""
        self.base_count = 1000

        self.sim.comp1.SA.Count = self.base_count
        self.sim.comp1.SB.Count = self.base_count
        self.sim.comp2.SC.Count = self.base_count

    def asserts(self, comp):

        SA, SB, SC = comp.SA.Count, comp.SB.Count, comp.SC.Count
        print(f"Comp: {comp.name}, SA: {SA}, SB: {SB}, SC: {SC}")

        self.assertEqual(SA, SB)
        self.assertEqual(SA + SC, self.base_count)
        if comp.name == "comp1":
            self.assertLessEqual(SA, 100)
            self.assertGreaterEqual(SA, 20)
        else:
            self.assertLessEqual(SA, 10)
            self.assertGreaterEqual(SA, 1)

    def run_sim(self):
        """Test after ENDTIME that the counts match with the results of a previous run taken as reference"""

        self.setUpSimulation()
        self.sim.newRun()
        self.setUpInitialConditions()

        self.sim.run(self.END_TIME)

        self.asserts(self.sim.comp1)
        self.asserts(self.sim.comp2)

    def test_tetopsplit_multi_comp_rd_n1(self):
        self.setUpMeshes("3tets_2patches_2comp.msh")
        self.run_sim()

    def test_tetopsplit_multi_comp_rd_n2(self):
        self.setUpMeshes("3tets_2patches_2comp.msh")
        self.run_sim()

    def test_tetopsplit_multi_comp_rd_split_mesh_n2(self):
        self.setUpMeshes("3tets_2patches_2comp_split2/3tets_2patches_2comp")
        self.run_sim()


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DistTetopsplitMultiCompRd, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
