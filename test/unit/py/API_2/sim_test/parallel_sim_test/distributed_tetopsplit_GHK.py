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

""" Unit tests for GHK reactions"""

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

FARADAY_CONSTANT = spc.physical_constants["Faraday constant"][0]
GAS_CONSTANT = spc.physical_constants["molar gas constant"][0]


class DistTetopsplitGHK(unittest.TestCase):
    """GHK channel unit test"""

    def setUp(self):
        self.setConstants()
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()

    def setConstants(self):
        # Equilibrium potential
        self.eqV = -65.0e-3

        self.DT = 1.0e-8

    def acr(self):
        """Concentrations and expected asymptotic concentration ratio"""
        acr = np.exp(-self.eqV * FARADAY_CONSTANT * self.SNa.valence / (GAS_CONSTANT * self.sim.Temp))

        return acr

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        self.model = Model()
        r = ReactionManager()
        with self.model:
            SNa = Species.Create(valence=2)

            state_0 = SubUnitState.Create()
            NaChan = Channel.Create([state_0])

            ssys = SurfaceSystem.Create()

            with ssys:
                NaCurr = GHKCurr.Create(NaChan[state_0], SNa, 1.0e-14, computeflux=True)

            self.ssys, self.NaCurr, self.SNa = ssys, NaCurr, SNa

    def setUpMeshes(self):
        """Create mesh and assign to it the VolumeSystem"""
        # "3tets_2patches_2comp_split2/3tets_2patches_2comp"
        self.mesh = DistMesh(os.path.join(MESH_DIR, "2_tets.msh"), scale=1e-6)
        with self.mesh:
            comp_i = Compartment.Create(conductivity=1.0)
            comp_o = Compartment.Create(conductivity=1.0)
            patch = Patch.Create(comp_i, comp_o, self.ssys)
            memb = Membrane.Create([patch], capacitance=0.0)

    def setUpSimulation(self):
        """Instantiate main simulator object"""
        rng = RNG('mt19937', 512, 7233)
        self.sim = Simulation('DistTetOpSplit', self.model, self.mesh, rng, isEfield=False)

    def setUpInitialConditions(self):
        """Set initial conditions"""
        state_0 = self.model.state_0
        self.sim.comp_i.SNa.Count = 20e6
        self.sim.comp_o.SNa.Count = 40e6
        self.sim.patch.NaChan[state_0].Count = 1

        self.conc_i_init = self.sim.comp_i.SNa.Conc * 1e3
        self.conc_o_init = self.sim.comp_o.SNa.Conc * 1e3

        self.sim.Temp = 273 + 30.0

    def test_tetopsplit_GHK_init_n1(self):

        """Test after ENDTIME that the counts match with the results of a previous run taken as reference"""

        NTESTS = 1000

        self.setUpInitialConditions()

        res = np.zeros((NTESTS))
        for i in range(NTESTS):
            print(f"test n: {i}/{NTESTS}")
            self.sim.newRun()
            self.sim.ALL(Membrane).Potential = self.eqV
            self.setUpInitialConditions()
            self.sim.run(self.DT)
            res[i] = self.sim.TRIS(self.mesh.patch.tris).NaCurr.I

        acr = self.acr()
        expected_init_current = (
            self.NaCurr.P
            * self.SNa.valence
            * FARADAY_CONSTANT
            * (-np.log(acr))
            * (self.conc_i_init - self.conc_o_init * acr)
            / (1.0 - acr)
        )

        self.assertTrue(np.isclose(expected_init_current, res.mean(), rtol=0.01, atol=0))

    def test_tetopsplit_GHK_asymptotic_n1(self):

        END_TIME = 10000 * self.DT
        self.setUpInitialConditions()
        self.sim.newRun()
        self.sim.ALL(Membrane).Potential = self.eqV
        self.setUpInitialConditions()
        self.sim.run(END_TIME)

        conc_i_asym = self.sim.comp_i.SNa.Conc * 1e3
        conc_o_asym = self.sim.comp_o.SNa.Conc * 1e3
        acr = self.acr()

        self.assertTrue(np.isclose(conc_i_asym, conc_o_asym * acr, rtol=0.01, atol=0))


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(DistTetopsplitGHK, "test"))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
