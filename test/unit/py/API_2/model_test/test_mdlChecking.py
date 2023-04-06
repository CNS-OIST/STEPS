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

""" Unit tests for model checking."""

import random
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.utils import *
import steps.simcheck as simcheck


class ComplexChecking(unittest.TestCase):
    """Test model checking."""
    def setUp(self, statesAsSpecies=True, order=NoOrdering):
        self.statesAsSpecies = statesAsSpecies
        self.order = order

        self.mdl = Model()
        self.r = ReactionManager()
        with self.mdl:
            S1, S2 = Species.Create()

            A1, A2, A3, B1, B2, B3 = SubUnitState.Create()

            SA, B = SubUnit.Create(
                [A1, A2, A3],
                [B1, B2, B3]
            )

            CC = Complex.Create(
                [SA, SA, SA, B, B], statesAsSpecies=self.statesAsSpecies, order=self.order
            )

            vsys1, vsys2 = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

        self.geom = Geometry()
        with self.geom:
            comp1, comp2 = Compartment.Create(
                Params(vsys1, 1), 
                Params(vsys2, 1)
            )
            patch = Patch.Create(comp1, comp2, ssys)

        self.S1, self.S2, self.A1, self.A2, self.A3, self.B1, self.B2, self.B3 = \
            S1, S2, A1, A2, A3, B1, B2, B3
        self.C = CC
        self.vsys1, self.vsys2, self.ssys = vsys1, vsys2, ssys
        self.comp1, self.comp2, self.patch = comp1, comp2, patch

        rng = RNG('mt19937', 512, 1234)

        self.sim = Simulation('Wmdirect', self.mdl, self.geom, rng)

    def testInvalidComplexReactions(self):
        S1, S2, A1, A2, A3, B1, B2, B3 =\
            self.S1, self.S2, self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        C = self.C
        r = self.r

        with self.mdl:
            with self.vsys1:
                C[A1, A1, A1, B1, B1] >r[1]> C[A3, A3, A3, B3, B3]
                r[1].K = 1
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertFalse(ok)
            self.assertEqual(len(errors), 1)
            self.assertEqual(len(warnings), 0)

            with self.vsys1:
                C[A3, A3, A3, B3, B3] + C[A1, A1, A1, B1, B1] >r[1]> C[A1, A1, A1, B1, B1] + C[A3, A3, A3, B3, B3]
                r[1].K = 10
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertFalse(ok)
            self.assertEqual(len(errors), 1)
            self.assertEqual(len(warnings), 0)

            with self.vsys1:
                C[A3, A3, A3, B3, B3] >r[1]> C[A1, A1, A1, B1, B1]
                r[1].K = 15
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertTrue(ok)
            self.assertEqual(len(errors), 0)
            self.assertEqual(len(warnings), 0)

            with self.vsys1:
                S1 >r[1]> S2
                r[1].K = 5
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertTrue(ok)
            self.assertEqual(len(errors), 0)
            self.assertEqual(len(warnings), 1)

            with self.vsys1:
                S2 >r[1]> S1
                r[1].K = 7
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertTrue(ok)
            self.assertEqual(len(errors), 0)
            self.assertEqual(len(warnings), 0)

            with self.ssys:
                C[A1, A1, A1, B1, B1].i >r[1]> C[A1, A1, A1, B1, B1].s
                r[1].K = 8
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertFalse(ok)
            self.assertEqual(len(errors), 1)
            self.assertEqual(len(warnings), 0)

            with self.ssys:
                C[A1, A1, A1, B1, B1].s >r[1]> C[A1, A1, A1, B1, B1].i
                r[1].K = 6
            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertTrue(ok)
            self.assertEqual(len(errors), 0)
            self.assertEqual(len(warnings), 0)

            with self.vsys1:
                for i in range(simcheck.RATE_NBREACS_THRESH):
                    S2 >r[1]> S1
                    r[1].K = random.randint(5, 10)

                S2 >r[1]> S1
                r[1].K = 100000000

                S2 >r[1]> S1
                r[1].K = 0

            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertTrue(ok)
            self.assertEqual(len(errors), 0)
            self.assertEqual(len(warnings), 2)

            with self.vsys1, self.ssys:
                (S1.i >r[4]> S1.o) + S2.o >r[5]> 2*S2
                r[4].K = 1
                r[5].K = 1

            ok, errors, warnings = simcheck.Check(self.sim, printmsgs = False)
            self.assertFalse(ok)
            self.assertEqual(len(errors), 1)
            self.assertEqual(len(warnings), 2)



def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ComplexChecking, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
