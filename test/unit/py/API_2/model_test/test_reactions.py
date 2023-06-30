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

import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.utils import *
import steps.simcheck


class ReactionDeclaration(unittest.TestCase):
    """Test reaction declaration."""
    def setUp(self):
        self.mdl = Model()
        self.r = ReactionManager()
        with self.mdl:
            S1, S2 = Species.Create()

            # A1, A2, A3, B1, B2, B3 = SubUnitState.Create()

            # SA, B = SubUnit.Create(
                # [A1, A2, A3],
                # [B1, B2, B3]
            # )

            # CC = Complex.Create(
                # [SA, SA, SA, B, B], statesAsSpecies=self.statesAsSpecies, order=self.order
            # )

            vsys1, vsys2 = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

        self.geom = Geometry()
        with self.geom:
            comp1, comp2 = Compartment.Create(
                Params(vsys1, 1), 
                Params(vsys2, 1)
            )
            patch = Patch.Create(comp1, comp2, ssys)

        self.S1, self.S2 = S1, S2
        self.vsys1, self.vsys2, self.ssys = vsys1, vsys2, ssys
        self.comp1, self.comp2, self.patch = comp1, comp2, patch

    def testReactionDeclaration(self):
        mdl, vsys1, vsys2, ssys, S1, S2, r = self.mdl, self.vsys1, self.vsys2, self.ssys, self.S1, self.S2, self.r

        with mdl:
            with vsys1:
                reac = Reaction.Create()
                S1 >reac> S2
                reac.K = 1
                with self.assertRaises(Exception):
                    S1 >reac> 2*S2

                with self.assertRaises(Exception):
                    S1 >r[13]< S2

                with self.assertRaises(Exception):
                    S1 >r[14]> 42

                with self.assertRaises(Exception):
                    S1 >r[14]> 'test'

                with self.assertRaises(Exception):
                    S1 >r[15]> S2
                    r[15].K = 'test'

                with self.assertRaises(Exception):
                    S1.o + S2.s >r[1]> S2.s + S1.i
                    r[1].K = 1

                S1 >r[18]
                with self.assertRaises(Exception):
                    r[18].K = 1

                with self.assertRaises(Exception):
                    r[19].K = 1

                S1 >r[20]> S2
                with self.assertRaises(Exception):
                    r[20].K = 1, 1

                S1 <r[21]> S2
                with self.assertRaises(Exception):
                    r[21].K = 1

                S1 <r[23]> S2
                reac = r[23]
                r[23].K = 1,2

                self.assertEqual(reac['fwd'].K, 1)
                self.assertEqual(reac['bkw'].K, 2)

                reac['fwd']
                reac['bkw']
                with self.assertRaises(KeyError):
                    reac['test']

                reac['fwd'].K = 10
                reac['bkw'].K = 20

                self.assertEqual(reac['fwd'].K, 10)
                self.assertEqual(reac['bkw'].K, 20)
                self.assertEqual(reac.K, (10, 20))

                # Tests on reactants and reaction side:

                with self.assertRaises(TypeError):
                    S1 + 'S2'

                with self.assertRaises(TypeError):
                    'test' * S1

                with self.assertRaises(TypeError):
                    -2 * S1
                
                with self.assertRaises(Exception):
                    None + (S1 + S2)

                with self.assertRaises(TypeError):
                    S1 + S2 + 'S2'

                with self.assertRaises(SyntaxError):
                    (S1 + S2) | S1



            with vsys1:
                S1 >r[22]> S2
                reac = r[22]
                r[22].K = 1
            with self.assertRaises(Exception):
                reac.K = 1

            with vsys1, ssys:
                (2 * S1.o >r[2]> S1.i) + S2.i >r[3]> 2*S2
                r[2].K = 1
                r[3].K = 1

                # Setting the location of several species at once
                Surf(S1 + S2)

                with self.assertRaises(Exception):
                    (S1.o >r[4]> S1.s) + S2.i >r[5]> 2*S2
                    r[4].K = 1
                    r[5].K = 1

                with self.assertRaises(Exception):
                    (S1.o >r[6]> S1.s) + S2.s >r[7]> 2*S2
                    r[6].K = 1
                    r[7].K = 1

            with self.assertRaises(Exception):
                with vsys1, ssys:
                    S1 + S2.o <r[45]> 2 * S1.i

            with ssys:
                with self.assertRaises(Exception):
                    S1 + S2 >r[11]> 2*S2
                    r[11].K = 1

            with self.assertRaises(TypeError):
                Surf('test')


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ReactionDeclaration, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

