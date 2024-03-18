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

""" Test to verify that we correctly call `create_dist_for_variable_sized` in diffusions.hpp when different ranks have
different number of species. Before we were checking only that each rank, locally, had every element with the same
number of species. We also need to check that that number is the same across all ranks. This test hangs if at least one
process does not call `create_dist_for_variable_sized`"""

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


class DistTetopsplitOmegahCreateDistForVariableSized(unittest.TestCase):
    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        self.mdl = Model()
        r = ReactionManager()

        with self.mdl:
            SA, SB, SC = Species.Create()
            vsys1, vsys2 = VolumeSystem.Create()
            with vsys1:
                SA > r[1] > SB
                r[1].K = 1
            with vsys2:
                SA > r[1] > SB + SC
                r[1].K = 1

            self.vsys1, self.vsys2 = vsys1, vsys2

    def setUpMeshes(self, path="two_comp_cyl.msh"):
        """Create mesh and assign to it the VolumeSystem"""
        scale = 1e-6
        path = os.path.join(MESH_DIR, path)

        self.mesh = DistMesh(path, scale=scale)

        with self.mesh:
            tets = TetList(tet for tet in self.mesh.tets if tet.center.z > 0)
            comp1 = Compartment.Create(tets, self.vsys1)
            comp2 = Compartment.Create(self.mesh.tets - tets, self.vsys2)
            self.comp1, self.comp2 = comp1, comp2

    def setUpSimulation(self, seed=7322):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, seed)
        self.sim = Simulation('DistTetOpSplit', self.mdl, self.mesh, rng)

    def test_tetopsplit_omegah_create_dist_for_variable_sized_n1(self):
        """ " This test should work in any case"""
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()

    def test_tetopsplit_omegah_create_dist_for_variable_sized_n2(self):
        """This test should hang in case a rank does not call `create_dist_for_variable_sized` in diffusions.hpp"""
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DistTetopsplitOmegahCreateDistForVariableSized, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
