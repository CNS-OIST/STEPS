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

import os
import unittest

import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class DistTetopsplitSRunit(unittest.TestCase):
    """Test simple reaction/diffusion system with DistTetopsplit"""

    def setUp(self):
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        self.model = Model()
        r = ReactionManager()
        with self.model:
            SA, SB, SC, SD = Species.Create()
            vsysi = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()
            with vsysi:
                SA + SB > r[1] > None
                r[1].K = 1e8

            with ssys:
                3 * SD.s  + 2 * SC.o > r[2] > 4 * SB.i
                r[2].K = 1e21

            self.vsysi, self.ssys = vsysi, ssys

    def setUpMeshes(self):
        """ Create mesh and assign to it the VolumeSystem"""
        self.mesh = DistMesh(os.path.join(MESH_DIR, "2_tets.msh"), scale=1e-6)
        with self.mesh:
            comp_i = Compartment.Create(vsys=self.vsysi)
            comp_o = Compartment.Create()
            patch = Patch.Create(comp_i, comp_o, self.ssys)

    def setUpSimulation(self):
        """ Instantiate main simulator object"""
        rng = RNG('mt19937', 512, 7233)
        self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng)

    def setUpInitialConditions(self):
        """ Set initial conditions """
        self.sim.comp_i.SA.Count = 4000
        self.sim.comp_i.SB.Count = 0
        self.sim.comp_o.SC.Count = 2000
        self.sim.patch.SD.Count = 6000

    def test_tetopsplit_sr_unit_n1(self):
        """ Test after ENDTIME that the counts match with the results of a previous run taken as reference """

        ENDTIME = 20
        expected_counts = [0.0, 0.0, 0.0, 3000.0]

        self.sim.newRun()
        self.setUpInitialConditions()
        self.sim.run(ENDTIME)

        counts = [
            self.sim.comp_i.SA.Count,
            self.sim.comp_i.SB.Count,
            self.sim.comp_o.SC.Count,
            self.sim.patch.SD.Count,
        ]
        print(f"\nFinal counts:    {''.join([f'{i:<8}' for i in counts])}")
        print(f"Expected counts: {''.join([f'{i:<8}' for i in expected_counts])}")
        self.assertEqual(counts, expected_counts)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistTetopsplitSRunit))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
