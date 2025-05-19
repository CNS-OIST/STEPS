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

""" Unit tests for diffusions only"""

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


class DistTetopsplitDiffOnly(unittest.TestCase):
    """Test only diffusions"""

    def setUp(self):
        self.setConstants()
        self.setUpModel()

    def setConstants(self):
        self.END_TIME = 0.001
        self.NTESTS = 1000

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        self.model = Model()
        r = ReactionManager()

        with self.model:
            vsys1 = VolumeSystem.Create()
            SA = Species.Create()
            with vsys1:
                DA1 = Diffusion.Create(SA, 1e-8)

        self.vsys1 = vsys1

    def setUpMeshes(self, path="diamond.msh"):
        """Create mesh and assign to it the VolumeSystem"""
        scale = 1e-6
        path = os.path.join(MESH_DIR, path)

        if self.steps_version == 4:
            self.mesh = DistMesh(path, scale=scale)
            with self.mesh:
                comp1 = Compartment.Create(vsys=self.vsys1)
        else:
            self.mesh = TetMesh.LoadGmsh(path, scale=scale)
            with self.mesh:
                comp1 = Compartment.Create(self.mesh.tets, self.vsys1)

    def setUpSimulation(self, seed=7322):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, seed)

        if self.steps_version == 4:
            self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng, isEfield=False)
        else:
            part = LinearMeshPartition(self.mesh, MPI.nhosts, 1, 1)
            self.sim = Simulation('TetOpSplit', self.model, self.mesh, rng, MPI.EF_NONE, part)

    def setUpInitialConditions(self):
        """Set initial conditions"""
        self.base_count = 1000
        self.injtetID = 0

        self.sim.TET(self.injtetID).SA.Count = self.base_count

    def run_sims(self):
        self.setUpMeshes()
        self.setUpSimulation()
        rs = ResultSelector(self.sim)
        counts = rs.TETS().SA.Count
        self.sim.toSave(counts)
        for i in range(self.NTESTS):
            self.sim.newRun()
            self.setUpInitialConditions()
            self.sim.run(self.END_TIME)
            counts.save()
            self.assertEqual(np.sum(counts.data[-1]), self.base_count)
        return counts.data

    def test_tetopsplit_diff_only_n1(self):
        """Test after ENDTIME that STEPS3 and 4 match over many tests"""

        self.steps_version = 3
        res3 = np.array(self.run_sims())
        self.assertTrue(np.all(res3 >= 0))
        res3 = res3.mean(axis=0)

        self.steps_version = 4
        res4 = np.array(self.run_sims())
        self.assertTrue(np.all(res4 >= 0))
        res4 = res4.mean(axis=0)

        print(f"Avg. counts: \nSTEPS3: {res3}\nSTEPS4: {res4}")

        self.assertTrue(np.allclose(res3, res4, rtol=0.01, atol=2))


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistTetopsplitDiffOnly))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
