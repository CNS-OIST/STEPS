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

""" Unit tests for simple reaction/diffusions"""

import logging
import os
import platform
import unittest

import numpy

import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class DistTetopsplitBox(unittest.TestCase):
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
            SA, SB, SC, SD, SE, SF, SG, SH, SI, SJ = Species.Create()
            self.vsys = VolumeSystem(name="vsys")
            with self.vsys:
                SA + SB < r[1] > SC
                SC + SD < r[2] > SE
                SF + SG < r[3] > SH
                SH + SI < r[4] > SJ
                r[1].K = 1000e6, 100
                r[2].K = 100e6, 10
                r[3].K = 10e6, 1
                r[4].K = 1e6, 1

                # The diffusion rules
                D1 = Diffusion.Create(SA, 100e-12)
                D2 = Diffusion.Create(SB, 90e-12)
                D3 = Diffusion.Create(SC, 80e-12)
                D4 = Diffusion.Create(SD, 70e-12)
                D5 = Diffusion.Create(SE, 60e-12)
                D6 = Diffusion.Create(SF, 50e-12)
                D7 = Diffusion.Create(SG, 40e-12)
                D8 = Diffusion.Create(SH, 30e-12)
                D9 = Diffusion.Create(SI, 20e-12)
                D10 = Diffusion.Create(SJ, 10e-12)

    def setUpMeshes(self):
        """Create mesh and assign to it the VolumeSystem"""
        self.mesh = DistMesh(os.path.join(MESH_DIR, "cube.msh"), scale=1e-6)
        with self.mesh:
            comp1 = Compartment.Create(self.vsys)

    def setUpSimulation(self):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, 7233)
        self.sim = Simulation(
            'DistTetOpSplit', self.model, self.mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK
        )

    def setUpInitialConditions(self):
        """Set initial conditions"""
        MOLECULE_RATIO = 1
        N0A = 1000
        N0B = 2000
        N0C = 3000
        N0D = 4000
        N0E = 5000
        N0F = 6000
        N0G = 7000
        N0H = 8000
        N0I = 9000
        N0J = 10000
        self.sim.comp1.SA.Count = N0A * MOLECULE_RATIO
        self.sim.comp1.SB.Count = N0B * MOLECULE_RATIO
        self.sim.comp1.SC.Count = N0C * MOLECULE_RATIO
        self.sim.comp1.SD.Count = N0D * MOLECULE_RATIO
        self.sim.comp1.SE.Count = N0E * MOLECULE_RATIO
        self.sim.comp1.SF.Count = N0F * MOLECULE_RATIO
        self.sim.comp1.SG.Count = N0G * MOLECULE_RATIO
        self.sim.comp1.SH.Count = N0H * MOLECULE_RATIO
        self.sim.comp1.SI.Count = N0I * MOLECULE_RATIO
        self.sim.comp1.SJ.Count = N0J * MOLECULE_RATIO

    def testTetopsplitDistBox(self):
        """Test after ENDTIME that the counts match with the results of a previous run taken as reference"""
        ENDTIME = 0.5
        atol = 10  # we accept results +- 10 or with only 1% relative error
        rtol = 0.01
        if platform.system() == "Darwin":
            expected_counts = [34.0, 1034.0, 673.0, 707.0, 8293.0, 281.0, 1281.0, 6322.0, 1603.0, 17397.0]
        else:
            expected_counts = [42.0, 1042.0, 669.0, 711.0, 8289.0, 302.0, 1302.0, 6369.0, 1671.0, 17329.0]

        self.sim.newRun()

        self.setUpInitialConditions()
        self.sim.run(ENDTIME)

        counts = self.sim.comp1.ALL(Species).Count
        if MPI.rank == 0:
            print(f"\nFinal counts:    {counts}")
            print(f"Expected counts: {expected_counts}")
            self.assertTrue(numpy.allclose(counts, expected_counts, rtol=rtol, atol=atol))


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DistTetopsplitBox, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
