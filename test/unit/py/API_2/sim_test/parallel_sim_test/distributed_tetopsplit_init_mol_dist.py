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

""" Unit tests initial molecule distributions"""

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


class DistTetopsplitInitMolDist(unittest.TestCase):
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
            SC = Species.Create()
            with vsys:
                # we need a dummy reaction or diffusion to make STEPS put the molecules
                DC = Diffusion.Create(SC, 0)

        self.vsys = vsys

    def setUpMeshes(self, path="diamond.msh"):
        """Create mesh and assign to it the VolumeSystem"""
        scale = 1e-6
        path = os.path.join(MESH_DIR, path)

        self.mesh = DistMesh(path, scale=scale)
        with self.mesh:
            comp1 = Compartment.Create(vsys=self.vsys)

    def setUpSimulation(self, seed=7322):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, seed)
        self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng, isEfield=False)

    def setUpInitialConditions(self, distributionMethod):
        """Set initial conditions"""

        self.base_count = Params(1e5, distributionMethod)
        self.sim.comp1.SC.Count = self.base_count

    def test_tetopsplit_init_mol_multinomial_dist_n2(self):
        """Check the multinomial distribution"""

        self.sim.newRun()
        self.setUpInitialConditions(DistributionMethod.MULTINOMIAL)

        counts = np.array(self.sim.TETS().SC.Count)
        vols = np.array([tet.Vol for tet in self.mesh.tets])
        counts_sum = np.sum(counts)
        vols_sum = np.sum(vols)
        if MPI.rank == 0:
            print(f"Counts: {counts} (sum: {counts_sum} == {self.base_count.args[0]})")
            self.assertEqual(counts_sum, self.base_count.args[0])

            expected_counts = self.base_count.args[0] * vols / vols_sum
            chi_frac = np.sum(np.divide(np.square(counts - expected_counts), expected_counts))

            threshold = 18.48  # conf. lvl.: 99%
            print(f"Total chi ratio: {chi_frac} (< {threshold})")
            self.assertLessEqual(chi_frac, threshold)

    def test_tetopsplit_init_mol_uniform_dist_n2(self):
        """Check the uniform distribution"""

        self.sim.newRun()
        self.setUpInitialConditions(DistributionMethod.UNIFORM)

        counts = np.array(self.sim.TETS().SC.Count)
        vols = np.array([tet.Vol for tet in self.mesh.tets])
        counts_sum = np.sum(counts)
        vols_sum = np.sum(vols)
        if MPI.rank == 0:
            expected_counts = self.base_count.args[0] * vols / vols_sum
            print(
                f"Counts:          {counts} (sum: {counts_sum} == {self.base_count.args[0]})\nExpected counts: {expected_counts}"
            )
            self.assertEqual(counts_sum, self.base_count.args[0])
            self.assertTrue(np.allclose(expected_counts, counts, rtol=1e-3, atol=1))


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DistTetopsplitInitMolDist, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
