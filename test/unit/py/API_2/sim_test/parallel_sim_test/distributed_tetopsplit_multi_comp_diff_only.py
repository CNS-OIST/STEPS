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

""" Unit tests for multi compartment diffusions only"""

import steps.interface

import os
import unittest

import numpy as np
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class DistTetopsplitMultiCompDiffOnly(unittest.TestCase):
    """Test only diffusions"""

    def setUp(self):
        self.setConstants()
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()

    def setConstants(self):
        self.DT = 1e-5
        self.NSTEPS = 101
        self.END_TIME = self.DT * self.NSTEPS

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        self.model = Model()
        r = ReactionManager()

        with self.model:
            vsys = VolumeSystem.Create()
            SC, SD, SE = Species.Create()
            with vsys:
                DC = Diffusion.Create(SC, 1e-8)

        self.vsys = vsys

    def setUpMeshes(self, path="3_tets.msh"):
        """Create mesh and assign to it the VolumeSystem"""
        scale = 1e-6
        path = os.path.join(MESH_DIR, path)

        self.mesh = DistMesh(path, scale=scale)
        with self.mesh:
            comp_i = Compartment.Create(vsys=self.vsys)
            comp_o = Compartment.Create(vsys=self.vsys)
            diffb = DiffBoundary.Create(comp_i.surface & comp_o.surface)

    def setUpSimulation(self, seed=7322):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, seed)
        self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng, isEfield=False)

    def setUpInitialConditions(self):
        """Set initial conditions"""
        self.base_count = 0.3e6
        self.sim.comp_i.SC.Count = self.base_count
        self.sim.comp_o.SC.Count = self.base_count

    def test_tetopsplit_multi_comp_diff_only_n2(self):
        """Test that asymptotically molecules distribute as 1/3 / 2/3"""

        self.sim.newRun()
        self.setUpInitialConditions()

        Ci, Co = self.sim.comp_i.SC.Count, self.sim.comp_o.SC.Count
        if MPI.rank == 0:
            print("Initial counts:")
            print(f"comp_i count: {Ci} (Exp: {self.base_count})")
            print(f"comp_o count: {Co} (Exp: {self.base_count})")
        self.assertEqual(Ci, self.base_count)
        self.assertEqual(Co, self.base_count)

        self.sim.run(self.DT)
        Ci, Co = self.sim.comp_i.SC.Count, self.sim.comp_o.SC.Count
        if MPI.rank == 0:
            print(f"Counts after 1 DT at {self.DT}s (diffusion off):")
            print(f"comp_i count: {Ci} (Exp: {self.base_count})")
            print(f"comp_o count: {Co} (Exp: {self.base_count})")
        self.assertEqual(Ci, self.base_count)
        self.assertEqual(Co, self.base_count)

        self.sim.diffb.SC.DiffusionActive = True
        self.sim.run(self.END_TIME)
        Ci, Co = self.sim.comp_i.SC.Count, self.sim.comp_o.SC.Count
        if MPI.rank == 0:
            print(f"Final counts at {self.END_TIME}s (diffusion on):")
            print(f"comp_i count: {Ci} (Exp: ~{2*self.base_count/3})")
            print(f"comp_o count: {Co} (Exp: ~{4*self.base_count/3})")
        self.assertEqual(
            Ci + Co, 2 * self.base_count, "Error: the diffusion boundary lost molecules during diffusion!"
        )
        self.assertTrue(
            np.allclose([Ci, Co], [2 * self.base_count / 3, 4 * self.base_count / 3], rtol=0.1, atol=5)
        )


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistTetopsplitMultiCompDiffOnly))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
