####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolver

class setCompVolTestCase(unittest.TestCase):
    """
    Test if the setCompVol method for wmdirect and wmrssa correctly updates 
    the Ccst of surface reactions.
    """
    def setUp(self):
        mdl = smodel.Model()

        self.v1 = 1e-20
        self.v2 = 2e-20
        self.a1 = 3e-14

        self.kreac = 200.0
        self.ksreac = 100.0

        S1 = smodel.Spec('S1', mdl)
        S2 = smodel.Spec('S2', mdl)
        S1S2 = smodel.Spec('S1S2', mdl)

        vsys = smodel.Volsys('vsys', mdl)
        ssys = smodel.Surfsys('ssys', mdl)

        smodel.Reac('reac', vsys, lhs=[S1, S2], rhs=[S2, S2], kcst=self.kreac)

        smodel.SReac('sreac', ssys, ilhs=[S1], slhs=[S2], srhs=[S1S2], kcst=self.ksreac)

        geom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', geom)
        comp1.setVol(self.v1)
        comp1.addVolsys('vsys')

        comp2 = sgeom.Comp('comp2', geom)
        comp2.setVol(self.v2)
        comp1.addVolsys('vsys')

        patch = sgeom.Patch('patch', geom, comp1, comp2)
        patch.addSurfsys('ssys')
        patch.setArea(self.a1)

        self.mdl, self.geom, self.rng = mdl, geom, srng.create('mt19937',512)
        self.rng.initialize(1234)

    def testSetCompVolWmdirect(self):
        self.sim = ssolver.Wmdirect(self.mdl, self.geom, self.rng)
        self.sim.reset()

        self._testSetCompVol()

    def testSetCompVolWmrssa(self):
        self.sim = ssolver.Wmrssa(self.mdl, self.geom, self.rng)
        self.sim.reset()

        self._testSetCompVol()

    def _testSetCompVol(self):
        Avogad = 6.02214076e23

        # Check that the volume and the Ccst values are as expected
        self.assertEqual(self.sim.getCompVol('comp1'), self.v1)
        # volume reac
        self.assertAlmostEqual(
            self.sim.getCompReacC('comp1', 'reac'),
            self.kreac / (1e3 * self.v1 * Avogad)
        )
        # surface reac
        self.assertAlmostEqual(
            self.sim.getPatchSReacC('patch', 'sreac'),
            self.ksreac / (1e3 * self.v1 * Avogad)
        )

        # Double the volume of comp1
        self.sim.setCompVol('comp1', 2 * self.v1)

        # Check that the volume was doubled
        self.assertEqual(self.sim.getCompVol('comp1'), 2 * self.v1)

        # Check that the Ccst was halved
        # volume reac
        self.assertAlmostEqual(
            self.sim.getCompReacC('comp1', 'reac'),
            self.kreac / (1e3 * 2 * self.v1 * Avogad)
        )
        # surface reac
        self.assertAlmostEqual(
            self.sim.getPatchSReacC('patch', 'sreac'),
            self.ksreac / (1e3 * 2 * self.v1 * Avogad)
        )

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(setCompVolTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
