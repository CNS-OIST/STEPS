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

""" Unit tests for DistTetOpSplit.setVertIClamp(...)"""

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


class DistTetopsplitSetVertIClamp(unittest.TestCase):
    """Test that vertices current injection works as expected with DistTetopsplit"""

    def setUp(self):
        self.Iinj = 1e-15
        self.capacitance = 0.01
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA = Species.Create()
            vsys = VolumeSystem.Create()
            with vsys:
                SA >r[1]> None
                r[1].K = 1
        self.mdl = mdl

    def setUpMeshes(self):
        """ Create mesh and assign to it the VolumeSystem"""
        self.mesh = DistMesh(os.path.join(MESH_DIR, '2_tets.msh'), scale=1e-6)
        with self.mesh:
            __MESH__ = Compartment.Create(self.mdl.vsys)
            __MESH_BOUNDARY__ = Patch.Create(__MESH__, None)
            membrane = Membrane.Create([__MESH_BOUNDARY__], capacitance=self.capacitance)
            __MESH__.Conductivity = 1

        self.area = __MESH_BOUNDARY__.Area
        self.nbVerts = len(self.mesh.verts)

    def setUpSimulation(self):
        """ Instantiate main simulator object"""
        rng = RNG('mt19937', 512, 7233)
        self.sim = Simulation('DistTetOpSplit', self.mdl, self.mesh, rng)

    def test_disttetopsplit_setVertIClamp_n2(self):
        """ Test after ENDTIME that the potential corresponds to what is expected given the current that should be injected"""
        ENDTIME = 0.1

        self.sim.newRun()

        self.sim.membrane.Potential = 0
        self.sim.VERTS().IClamp = self.Iinj

        self.sim.run(ENDTIME)

        avgPot = sum(self.sim.VERTS().V) / self.nbVerts
        expectedPot = ENDTIME * (self.Iinj * self.nbVerts) / (self.capacitance * self.area)

        if MPI.rank == 0:
            print(f'AvgPot = {avgPot*1e3}mV, expectedPot={expectedPot*1e3}mV, difference: {abs(avgPot- expectedPot)*1e3}mV')
            self.assertAlmostEqual(avgPot, expectedPot)

    def test_disttetopsplit_setVertIClamp_n1(self):
        # Check that we have the same behavior for 1 rank
        self.test_disttetopsplit_setVertIClamp_n2()


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DistTetopsplitSetVertIClamp, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
