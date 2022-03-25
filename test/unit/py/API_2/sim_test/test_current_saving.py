####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

""" Unit tests for current saving."""

import os
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

from . import base_model


class TetCurrentSaving(base_model.TetTestModelFramework):
    """Test checkpoint and restore methods for Wmdirect solver"""
    def setUp(self):
        super().setUp()

        self.efielddt = 1e-3
        self.deltaT = 1e-3 + 1e-10

        self.initChan1Cl = 10
        self.Chan1_G = 20e-15
        self.Chan1_rev = -77e-3
        self.Chan1_P = 2.5e-15

        self.newMdl = self.get_API2_Mdl()
        self.newGeom = self.get_API2_Geom(self.newMdl)

        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Tetexact', nmdl, ngeom, nrng, True)
        nsim.EfieldDT = self.efielddt
        return nsim

    def init_API2_sim(self, sim):
        super().init_API2_sim(sim)

        sim.ALL(Compartment, Patch).ALL(Species).Clamped = True

        sim.ALL(Patch).Chan1[self.newMdl.chancl].Count = 0
        sim.ALL(Patch).Chan1[self.newMdl.chanop].Count = self.initChan1Cl
        sim.ALL(Patch).Chan1[self.newMdl.chancl].Clamped = True
        sim.ALL(Patch).Chan1[self.newMdl.chanop].Clamped = True


    def testSmallCurrentRecording(self):

        rs = ResultSelector(self.newSim)

        curr = rs.TRIS(self.newGeom.patch.tris).Chan1_GHK_I.I

        self.newSim.toSave(curr, dt=self.deltaT)

        self.newSim.newRun()
        self.init_API2_sim(self.newSim)

        self.newSim.run(self.endTime)

        if MPI._shouldWrite:
            nbZeroes = sum(1 for v in curr.data[0,:,0] if v == 0)
            self.assertLess(nbZeroes, 2)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TetCurrentSaving, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

