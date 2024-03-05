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

""" Unit tests for complex reaction declaration."""

import unittest

from steps import interface

from steps.model import *

class ComplexCreation(unittest.TestCase):
    """Test Complex creation."""
    def setUp(self):
        self.mdl = Model()
        self.statesAsSpecies = False

    def testNormalCreation(self):
        with self.mdl:
            A1, A2, A3, B1, B2, B3 = SubUnitState.Create()

            SA, B = SubUnit.Create(
                [A1, A2, A3],
                [B1, B2, B3]
            )

            CC = Complex.Create(
                [SA, SA, SA, B, B], statesAsSpecies=self.statesAsSpecies, order=NoOrdering
            )

            with self.assertRaises(Exception):
                CC2 = Complex.Create(
                    [SA, SA, SA, B, B], statesAsSpecies=self.statesAsSpecies, order=StrongOrdering
                )

    def testCreactionFromSubStates(self):
        with self.mdl:
            C1, C2, C3 = SubUnitState.Create()

            CC = Complex.Create(
                [C1, C2, C3], 
                statesAsSpecies=self.statesAsSpecies, order=NoOrdering
            )

            with self.assertRaises(Exception):
                CC2 = Complex.Create(
                    [C1, C2, C3], 
                    statesAsSpecies=self.statesAsSpecies, order=StrongOrdering
                )

class ComplexCreationStatesAsSpecies(unittest.TestCase):
    """Test Complex creation."""
    def setUp(self):
        self.mdl = Model()
        self.statesAsSpecies = True
        self.order = NoOrdering

    def testNormalCreation(self):
        with self.mdl:
            A1, A2, A3, B1, B2, B3 = SubUnitState.Create()

            SA, B = SubUnit.Create(
                [A1, A2, A3],
                [B1, B2, B3]
            )

            CC = Complex.Create(
                [SA, SA, SA, B, B], statesAsSpecies=self.statesAsSpecies, order=self.order
            )

    def testCreactionFromSubStates(self):
        with self.mdl:
            C1, C2, C3 = SubUnitState.Create()

            CC = Complex.Create(
                [C1, C2, C3], 
                statesAsSpecies=self.statesAsSpecies, order=self.order
            )

            # TODO write test to check that exception is raised when trying to create steps complexes with some ordering

class ComplexReacTest5Subunits(unittest.TestCase):
    """Base class for complex reaction testing."""
    def setUp(self, statesAsSpecies=True, order=NoOrdering):
        self.mdl = Model()
        self.r = ReactionManager()
        with self.mdl:
            S1, S2 = Species.Create()

            A1, A2, A3, B1, B2, B3, A1p, A2p, A3p, B1p, B2p, B3p, E1, E2 = SubUnitState.Create()

            SA, B, SAp, Bp, E = SubUnit.Create(
                [A1, A2, A3],
                [B1, B2, B3],
                [A1p, A2p, A3p],
                [B1p, B2p, B3p],
                [E1, E2]
            )

            CC = Complex.Create(
                [SA, SA, SA, B, B], statesAsSpecies=statesAsSpecies, order=order
            )

            CD = Complex.Create([SAp, Bp], statesAsSpecies=statesAsSpecies)

            F = Complex.Create([E, E], statesAsSpecies=statesAsSpecies)

            vsys = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

            self.S1, self.S2 = S1, S2
            self.A1, self.A2, self.A3, self.B1, self.B2, self.B3 = A1, A2, A3, B1, B2, B3
            self.A1p, self.A2p, self.A3p, self.B1p, self.B2p, self.B3p = A1p, A2p, A3p, B1p, B2p, B3p
            self.E1, self.E2 = E1, E2
            self.A, self.B, self.E = SA, B, E
            self.C, self.D, self.F = CC, CD, F
            self.vsys, self.ssys = vsys, ssys

    def _checkReactions(self, expected, actual, checkRates=False):
        """
        Check that the reactions that were declared match the ones that were expected.
        expected should be a list of tuples. Each tuple should have 2 elements corresponding to 
        lhs and rhs. Each element is a list of ComplexStates and/or Species.
        actual should be a list of steps reaction objects.
        """
        count = lambda v: [(v.count(e), e) for e in v]
        expected = [e if len(e) > 2 else e + (None,) for e in expected]

        expectedSet = {}
        for lhs, rhs, rate in expected:
            lhsNames, rhsNames = [], []
            for e in lhs:
                lhsNames += [so.getID() for so in e._getStepsObjects()]
            for e in rhs:
                rhsNames += [so.getID() for so in e._getStepsObjects()]
            key = (frozenset(count(lhsNames)), frozenset(count(rhsNames)))
            if key not in expectedSet:
                expectedSet[key] = rate
            elif checkRates:
                expectedSet[key] += rate

        actualSet = {}
        for reac in actual:
            lhsNames = [s.getID() for s in reac.getLHS()]
            rhsNames = [s.getID() for s in reac.getRHS()]
            rate = reac.getKcst() if checkRates else None
            key = (frozenset(count(lhsNames)), frozenset(count(rhsNames)))
            if key not in actualSet:
                actualSet[key] = rate
            elif checkRates:
                actualSet[key] += rate

        self.assertDictEqual(expectedSet, actualSet)

    def _checkSurfReactions(self, expected, actual, checkRates=False):
        """
        Check that the surface reactions that were declared match the ones that were expected.
        expected should be a list of tuples. Each tuple should have 2 elements corresponding to 
        lhs and rhs. Each element is a list of ReactingElements for ComplexStates and/or Species.
        actual should be a list of steps reaction objects.
        """
        count = lambda v: [(v.count(e), e) for e in v]
        expected = [e if len(e) > 2 else e + (None,) for e in expected]

        expectedSet = set()
        for lhs, rhs, rate in expected:
            lhsNames, rhsNames = [], []
            for e in lhs:
                lhsNames += [(e.loc, so.getID()) for so in e._elem._getStepsObjects()]
            for e in rhs:
                rhsNames += [(e.loc, so.getID()) for so in e._elem._getStepsObjects()]
            if checkRates:
                expectedSet.add((frozenset(count(lhsNames)), frozenset(count(rhsNames)), rate))
            else:
                expectedSet.add((frozenset(count(lhsNames)), frozenset(count(rhsNames))))

        actualSet = set()
        for reac in actual:
            lhsNames =  [(Location.IN, s.getID()) for s in reac.getILHS()]
            rhsNames =  [(Location.IN, s.getID()) for s in reac.getIRHS()]
            lhsNames += [(Location.SURF, s.getID()) for s in reac.getSLHS()]
            rhsNames += [(Location.SURF, s.getID()) for s in reac.getSRHS()]
            lhsNames += [(Location.OUT, s.getID()) for s in reac.getOLHS()]
            rhsNames += [(Location.OUT, s.getID()) for s in reac.getORHS()]
            rate = reac.getKcst()
            if checkRates:
                actualSet.add((frozenset(count(lhsNames)), frozenset(count(rhsNames)), rate))
            else:
                actualSet.add((frozenset(count(lhsNames)), frozenset(count(rhsNames))))

        self.assertEqual(expectedSet, actualSet)


class CompStateGeneral(ComplexReacTest5Subunits):
    """Test complex indexing."""
    def setUp(self, statesAsSpecies=True):
        super().setUp(statesAsSpecies, order=StrongOrdering)

    def testIndexing(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C = self.C

        self.assertTrue(A1 in (A1|A2))
        self.assertFalse(A1 in (A2|A3))
        self.assertFalse(A1 in ~A1)

        cs = C[A1, A2, A3, B1, B2]
        self.assertIsInstance(cs, ComplexState)
        self.assertEqual(cs, C[A1, A2, A3, B1, B2])

        cs = C[A1, :, A3, B1, B2]
        #self.assertIsInstance(cs, ComplexSelector)
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A1, A3, B1, B2],
                C[A1, A2, A3, B1, B2],
                C[A1, A3, A3, B1, B2],
            ])
        )

        cs = C[A1, A2, ..., B2]
        #self.assertIsInstance(cs, ComplexSelector)
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A1, B1, B2],
                C[A1, A2, A2, B1, B2],
                C[A1, A2, A3, B1, B2],
                C[A1, A2, A1, B2, B2],
                C[A1, A2, A2, B2, B2],
                C[A1, A2, A3, B2, B2],
                C[A1, A2, A1, B3, B2],
                C[A1, A2, A2, B3, B2],
                C[A1, A2, A3, B3, B2],
            ])
        )

        cs = C[..., A1|A3, A3, B1, B2]
        #self.assertIsInstance(cs, ComplexSelector)
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A1, A3, B1, B2],
                C[A1, A3, A3, B1, B2],
                C[A2, A1, A3, B1, B2],
                C[A2, A3, A3, B1, B2],
                C[A3, A1, A3, B1, B2],
                C[A3, A3, A3, B1, B2],
            ])
        )

        C[...]

        with self.assertRaises(Exception):
            A1|B1

        with self.assertRaises(Exception):
            A1|42

        with self.assertRaises(Exception):
            C[B1, A2, A3, B1, B2]

        with self.assertRaises(Exception):
            C[A1, A2, A3, 42, B2]

        with self.assertRaises(Exception):
            C[A1, :, B3, B1, B2]

        with self.assertRaises(Exception):
            C[A1, A2, A3, B1, B2, B1]

        with self.assertRaises(Exception):
            C[~(A1|A2|A3), A2, A3, B1, B2]

        with self.assertRaises(Exception):
            C[..., A1, ..., B2]

        with self.assertRaises(Exception):
            C[A1, A2, A3]

        with self.assertRaises(Exception):
            A1[None]

        with self.assertRaises(Exception):
            A1[[]]

        with self.assertRaises(Exception):
            with C[:, A2, :, B1, B2] as C2:
                A1[1, C2, 2]

        with self.assertRaises(Exception):
            A1[1, 1]

        with self.assertRaises(Exception):
            with C[:, A2, :, B1, B2] as C2:
                A1[C2,C2]

        with self.assertRaises(Exception):
            with C[:, A2, :, B1, B2] as C2:
                A1[C2]|A2

        with self.assertRaises(Exception):
            with C[:, A2, :, B1, B2] as C2, C[:, A2, :, B1, B2] as C3:
                (A1[C3]|A2[C3])[C2]

    def testCombining(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A1p, A2p, A3p, B1p, B2p, B3p = self.A1p, self.A2p, self.A3p, self.B1p, self.B2p, self.B3p
        A, B = self.A, self.B
        C, D = self.C, self.D

        cs = C[A1, :, A3, B1, B2] | C[A1, A2, A3, :, B2]
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A1, A3, B1, B2],
                C[A1, A2, A3, B1, B2],
                C[A1, A3, A3, B1, B2],
                C[A1, A2, A3, B1, B2],
                C[A1, A2, A3, B2, B2],
                C[A1, A2, A3, B3, B2],
            ])
        )

        cs = C[A1, :, A3, B1, B2] & C[A1, A2, A3, :, B2]
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A3, B1, B2],
            ])
        )

        cs = C[A1, :, A3|A2, ~B1, B2] & C[A1, A2|A1, ~A2, :, ~B3]
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A3, B2, B2],
                C[A1, A2, A3, B3, B2],
                C[A1, A1, A3, B2, B2],
                C[A1, A1, A3, B3, B2],
            ])
        )

        with self.assertRaises(Exception):
            C[A1, :, A3, B1, B2] & A2

        with self.assertRaises(Exception):
            C[A1, :, A3, B1, B2] | A2

        with self.assertRaises(Exception):
            C[A1, :, A3, B1, B2] | D[A1p, B1p]

        with self.assertRaises(Exception):
            C[A1, :, A3, B1, B2] & D[A1p, :]

        with self.assertRaises(Exception):
            A1[1]|A2[2]

        with self.assertRaises(Exception):
            C[A1, A1, A3, B1, B2] | 42

    def testInjecting(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A1p, A2p, A3p, B1p, B2p, B3p = self.A1p, self.A2p, self.A3p, self.B1p, self.B2p, self.B3p
        C, D = self.C, self.D

        cs = C[A1, ..., B1, B2] << A2 << A3
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A3, B1, B2],
                C[A1, A3, A2, B1, B2],
            ])
        )

        cs = C[A1, ..., B2] << A2 << (A3|A1)
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A3, B1, B2],
                C[A1, A2, A1, B1, B2],
                C[A1, A2, A3, B2, B2],
                C[A1, A2, A1, B2, B2],
                C[A1, A2, A3, B3, B2],
                C[A1, A2, A1, B3, B2],
                C[A1, A3, A2, B1, B2],
                C[A1, A1, A2, B1, B2],
                C[A1, A3, A2, B2, B2],
                C[A1, A1, A2, B2, B2],
                C[A1, A3, A2, B3, B2],
                C[A1, A1, A2, B3, B2],
            ])
        )

        cs = C[A1, ..., B2] << A2[1] << A3 << (B1|B3)
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A3, B1, B2],
                C[A1, A2, A3, B3, B2],
            ])
        )

        cs = C[A1, ...] << A2[1] << A3 << (B1|B3)[3]
        self.assertEqual(
            set(cs), 
            set([
                C[A1, A2, A3, B1, B1],
                C[A1, A2, A3, B3, B1],
                C[A1, A2, A3, B1, B2],
                C[A1, A2, A3, B3, B2],
                C[A1, A2, A3, B1, B3],
                C[A1, A2, A3, B3, B3],
            ])
        )

        cs = C[..., B1, B1] << 2 * A1 + A2
        self.assertEqual(
            set(cs), 
            set([
                C[A2, A1, A1, B1, B1],
                C[A1, A2, A1, B1, B1],
                C[A1, A1, A2, B1, B1],
            ])
        )

        with self.assertRaises(Exception):
            C[...] << 42

        with self.assertRaises(Exception):
            C[A1, A2, A3, :, B2] << A1

        with self.assertRaises(Exception):
            C[...] << A1.s

        with self.assertRaises(Exception):
            C[...] << 4 * A1

    def testIdentifiers(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A1p, A2p, A3p, B1p, B2p, B3p = self.A1p, self.A2p, self.A3p, self.B1p, self.B2p, self.B3p
        C, D = self.C, self.D
        r = self.r

        with self.mdl:
            with self.vsys:
                with C[:, :, :, B1, B1]:
                    A1['a'] + A2['b'] >r[1]> A2['a'] + A3['b']
                r[1].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B1]], [C[A2, A3, A1, B1, B1]], 1),
                ([C[A1, A2, A2, B1, B1]], [C[A2, A3, A2, B1, B1]], 1),
                ([C[A1, A2, A3, B1, B1]], [C[A2, A3, A3, B1, B1]], 1),
                ([C[A1, A1, A2, B1, B1]], [C[A2, A1, A3, B1, B1]], 1),
                ([C[A1, A2, A2, B1, B1]], [C[A2, A2, A3, B1, B1]], 1),
                ([C[A1, A3, A2, B1, B1]], [C[A2, A3, A3, B1, B1]], 1),
                ([C[A1, A1, A2, B1, B1]], [C[A1, A2, A3, B1, B1]], 1),
                ([C[A2, A1, A2, B1, B1]], [C[A2, A2, A3, B1, B1]], 1),
                ([C[A3, A1, A2, B1, B1]], [C[A3, A2, A3, B1, B1]], 1),
                ([C[A2, A1, A1, B1, B1]], [C[A3, A2, A1, B1, B1]], 1),
                ([C[A2, A1, A2, B1, B1]], [C[A3, A2, A2, B1, B1]], 1),
                ([C[A2, A1, A3, B1, B1]], [C[A3, A2, A3, B1, B1]], 1),
                ([C[A2, A1, A1, B1, B1]], [C[A3, A1, A2, B1, B1]], 1),
                ([C[A2, A2, A1, B1, B1]], [C[A3, A2, A2, B1, B1]], 1),
                ([C[A2, A3, A1, B1, B1]], [C[A3, A3, A2, B1, B1]], 1),
                ([C[A1, A2, A1, B1, B1]], [C[A1, A3, A2, B1, B1]], 1),
                ([C[A2, A2, A1, B1, B1]], [C[A2, A3, A2, B1, B1]], 1),
                ([C[A3, A2, A1, B1, B1]], [C[A3, A3, A2, B1, B1]], 1),
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects(), checkRates=True)


class CompStateStrongOrdering(ComplexReacTest5Subunits):
    """Test equality of states under strong ordering."""
    def setUp(self, statesAsSpecies=True):
        super().setUp(statesAsSpecies, order=StrongOrdering)

    def testStateOrdering(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C = self.C
        self.assertEqual(C[A1, A2, A3, B1, B2], C[A1, A2, A3, B1, B2])
        self.assertNotEqual(C[A1, A2, A3, B1, B2], C[A2, A1, A3, B1, B2])


class CompStateNoOrdering(ComplexReacTest5Subunits):
    """Test equality of states under no ordering."""
    def setUp(self, statesAsSpecies=True):
        super().setUp(statesAsSpecies, order=NoOrdering)

    def testStateOrdering(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C = self.C
        self.assertEqual(C[A1, A2, A3, B1, B2], C[A1, A2, A3, B1, B2])
        self.assertEqual(C[A1, A2, A3, B1, B2], C[A3, A2, A1, B2, B1])
        self.assertEqual(C[A1, A1, A3, B1, B2], C[A1, A3, A1, B2, B1])
        self.assertNotEqual(C[A1, A2, A3, B1, B2], C[A2, A2, A3, B1, B2])
        self.assertNotEqual(C[A1, A2, A3, B1, B2], C[A1, A2, A3, B1, B1])


class CompStateRotationalSymmetryOrdering(ComplexReacTest5Subunits):
    """Test equality of states under rotational symetry ordering."""
    def setUp(self, statesAsSpecies=True):
        func = lambda s: RotationalSymmetryOrdering(s[:3]) + RotationalSymmetryOrdering(s[3:])
        super().setUp(statesAsSpecies, order=func)

    def testStateOrdering(self):
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C = self.C
        self.assertEqual(C[A1, A2, A3, B1, B2], C[A1, A2, A3, B1, B2])
        self.assertEqual(C[A1, A2, A3, B1, B2], C[A2, A3, A1, B2, B1])
        self.assertEqual(C[A1, A2, A3, B1, B2], C[A3, A1, A2, B2, B1])
        self.assertEqual(C[A1, A1, A3, B1, B2], C[A1, A3, A1, B2, B1])
        self.assertNotEqual(C[A1, A2, A3, B1, B2], C[A1, A3, A2, B1, B2])
        self.assertNotEqual(C[A1, A2, A3, B1, B2], C[A1, A2, A3, B1, B1])


class CompSelNoOrdReacTestCase(ComplexReacTest5Subunits):
    """Test reaction declaration with Complex Selectors and no ordering."""
    def testFullyDeterminedFwdReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A1p, A2p, A3p, B1p, B2p, B3p = self.A1p, self.A2p, self.A3p, self.B1p, self.B2p, self.B3p
        A, B, C, D = self.A, self.B, self.C, self.D
        C2, C3, D2 = self.C.get(), self.C.get(), self.D.get()
        with self.mdl:
            with self.vsys:
                C2[A1, A2, A3, B1, B2] >r[1]> C2[A2, A2, A3, B1, B2]
                r[1].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects())

            with self.vsys:
                C2[A1, A2, A3, B1, B2] | C2[A2, A2, A2, B3, B2] >r[2]> C2[A2, A2, A3, B1, B2]
                r[2].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B1, B2]]),
                ([C2[A2, A2, A2, B3, B2]], [C2[A2, A2, A3, B1, B2]]),
            ]
            self._checkReactions(expectedReacs, r[2]._getStepsObjects())

            with self.vsys:
                C2[A1, A2, A3, B1, B2] + C3[A1, A1, A2, B3, B3] >r[5]> C2[A2, A2, A3, B1, B2] + C3[A1, A1, A1, B1, B1]
                r[5].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2], C3[A1, A1, A2, B3, B3]], [C2[A2, A2, A3, B1, B2], C3[A1, A1, A1, B1, B1]])
            ]
            self._checkReactions(expectedReacs, r[5]._getStepsObjects())

            with self.vsys:
                C2[A1, A2, A3, B1, B2] + D2[A1p, B3p] >r[6]> C2[A2, A2, A3, B1, B2] + D2[A3p, B1p]
                r[6].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2], D2[A1p, B3p]], [C2[A2, A2, A3, B1, B2], D2[A3p, B1p]])
            ]
            self._checkReactions(expectedReacs, r[6]._getStepsObjects())

            with self.vsys:
                C2[A1, A2, A3, B1, B2] + D2[A1p, B3p] >r[7]> C2[A2, A2, A3, B1, B2]
                r[7].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2], D2[A1p, B3p]], [C2[A2, A2, A3, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[7]._getStepsObjects())

            with self.vsys:
                C2[A1, A2, A3, B1, B2] + C2[A1, A2, A3, B1, B2] >r[8]> C2[A1, A1, A3, B1, B2] + C2[A1, A1, A3, B1, B2]
                r[8].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2], C2[A1, A2, A3, B1, B2]], [C2[A1, A1, A3, B1, B2], C2[A1, A1, A3, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[8]._getStepsObjects())

            with self.vsys:
                C2[A1, A2, A3, B1, B2] >r[9]> C2[A2, A2, A3, B1, B2] + (D2[A1p, ~B2p] & D2[A1p, ~B3p])
                r[9].K = 1
            expectedReacs = [
                ([C[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B1, B2], D[A1p, B1p]])
            ]
            self._checkReactions(expectedReacs, r[9]._getStepsObjects())


            with self.ssys:
                # Surface reactions
                C2[A1, A2, A3, B1, B2].s >r[3]> C2[A2, A2, A3, B1, B2].s
                r[3].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2].s], [C2[A2, A2, A3, B1, B2].s])
            ]
            self._checkSurfReactions(expectedReacs, r[3]._getStepsObjects())

            with self.ssys:
                C2[A1, A2, A3, B1, B2].i >r[4]> C2[A2, A2, A3, B1, B2].o
                r[4].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2].i], [C2[A2, A2, A3, B1, B2].o])
            ]
            self._checkSurfReactions(expectedReacs, r[4]._getStepsObjects())

    def testFullyDeterminedBidirReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C2 = self.C.get()
        with self.mdl:
            with self.vsys:
                C2[A1, A2, A3, B1, B2] <r[1]> C2[A2, A2, A3, B1, B2]
                r[1].K = 1, 2
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B1, B2]]),
                ([C2[A2, A2, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects())

    def testPartiallyDeterminedFwdReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C2 = self.C.get()
        with self.mdl:
            with self.vsys:
                C2[A1, A2, :, B1, B2] >r[1]> C2[A2, A2, :, B1, B2]
                r[1].K = 1
            expectedReacs = [
                ([C2[A1, A2, A1, B1, B2]], [C2[A2, A2, A1, B1, B2]], 1),
                ([C2[A1, A2, A2, B1, B2]], [C2[A2, A2, A2, B1, B2]], 1),
                ([C2[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B1, B2]], 1)
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects(), checkRates=True)

            with self.vsys:
                C2[A1, A2, A3, ~B1, B2] >r[2]> C2[A2, A2, A3, ~(B2|B3), B2]
                r[2].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B2, B2]], [C2[A2, A2, A3, B1, B2]], 1),
                ([C2[A1, A2, A3, B3, B2]], [C2[A2, A2, A3, B1, B2]], 1),
            ]
            self._checkReactions(expectedReacs, r[2]._getStepsObjects(), checkRates=True)

            with self.vsys:
                C2[A1, A2, :, B1|B3, B2] >r[3]> C2[A2, A2, :, B2, B2]
                r[3].K = 1
            expectedReacs = [
                ([C2[A1, A2, A1, B1, B2]], [C2[A2, A2, A1, B2, B2]], 1),
                ([C2[A1, A2, A2, B1, B2]], [C2[A2, A2, A2, B2, B2]], 1),
                ([C2[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B2, B2]], 1),
                ([C2[A1, A2, A1, B3, B2]], [C2[A2, A2, A1, B2, B2]], 1),
                ([C2[A1, A2, A2, B3, B2]], [C2[A2, A2, A2, B2, B2]], 1),
                ([C2[A1, A2, A3, B3, B2]], [C2[A2, A2, A3, B2, B2]], 1)
            ]
            self._checkReactions(expectedReacs, r[3]._getStepsObjects(), checkRates=True)

            with self.vsys:
                C2[A1, A2, A3, B1, :] | C2[:, A1, A3, B1, B2] >r[6]> C2[A1, A2, A3, B1, B2]
                r[6].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B1]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A1, A2, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A1, A2, A3, B1, B3]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A1, A1, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A3, A1, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 1),
            ]
            self._checkReactions(expectedReacs, r[6]._getStepsObjects(), checkRates=True)

            with self.vsys:
                C2[A1, A1, A3, B1, B2] | C2[A1, A2, A3, B1, :]  >r[7]> C2[A1, A2, A3, B1, B2]
                r[7].K = 1
            expectedReacs = [
                ([C2[A1, A2, A3, B1, B1]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A1, A2, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A1, A2, A3, B1, B3]], [C2[A1, A2, A3, B1, B2]], 1),
                ([C2[A1, A1, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 1),
            ]
            self._checkReactions(expectedReacs, r[7]._getStepsObjects(), checkRates=True)

            with self.vsys:
                C2[~A1, :, :, B1, B2] >r[10]> C2[A1, :, :, B1, B2]
                r[10].K = 1
            expectedReacs = [
                ([C2[A2, A1, A1, B1, B2]], [C2[A1, A1, A1, B1, B2]], 1),
                ([C2[A2, A1, A2, B1, B2]], [C2[A1, A1, A2, B1, B2]], 1),
                ([C2[A2, A1, A3, B1, B2]], [C2[A1, A1, A3, B1, B2]], 0.5),
                ([C2[A2, A2, A2, B1, B2]], [C2[A1, A2, A2, B1, B2]], 1),
                ([C2[A2, A2, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 0.5),
                ([C2[A2, A3, A3, B1, B2]], [C2[A1, A3, A3, B1, B2]], 0.5),
                ([C2[A3, A1, A1, B1, B2]], [C2[A1, A1, A1, B1, B2]], 1),
                ([C2[A3, A1, A2, B1, B2]], [C2[A1, A1, A2, B1, B2]], 0.5),
                ([C2[A3, A1, A3, B1, B2]], [C2[A1, A1, A3, B1, B2]], 1),
                ([C2[A3, A2, A2, B1, B2]], [C2[A1, A2, A2, B1, B2]], 0.5),
                ([C2[A3, A2, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]], 0.5),
                ([C2[A3, A3, A3, B1, B2]], [C2[A1, A3, A3, B1, B2]], 1),
            ]
            self._checkReactions(expectedReacs, r[10]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with self.assertRaises(Exception):
                    C2[A1, A2, :, B1|B3, B2] >r[4]> C2[A2, A2, :, :, B2]

                C3 = self.C.get()
                with self.assertRaises(Exception):
                    C2[A1, A2, :, B1, B2] >r[5]> C3[A2, A2, :, B1, B2]

                with self.assertRaises(Exception):
                    C2[A1, A2, :, B1, B2] >r[8]> (C2[A2, A2, A1, B1, B2] | C2[A2, A2, A3, B1, B2])

                with self.assertRaises(Exception):
                    C2[A1, A2, :, B1, B2] >r[9]> C2[A2, A2, ~A1, B1, B2]


    def testPartiallyDeterminedBidirReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B = self.A, self.B
        C2 = self.C.get()
        with self.mdl:
            with self.vsys:
                C2[A1, A2, :, B1, B2] <r[1]> C2[A2, A2, :, B1, B2]
                r[1].K = 1, 1
            expectedReacs = [
                ([C2[A1, A2, A1, B1, B2]], [C2[A2, A2, A1, B1, B2]]),
                ([C2[A1, A2, A2, B1, B2]], [C2[A2, A2, A2, B1, B2]]),
                ([C2[A1, A2, A3, B1, B2]], [C2[A2, A2, A3, B1, B2]]),
                ([C2[A2, A2, A1, B1, B2]], [C2[A1, A2, A1, B1, B2]]),
                ([C2[A2, A2, A2, B1, B2]], [C2[A1, A2, A2, B1, B2]]),
                ([C2[A2, A2, A3, B1, B2]], [C2[A1, A2, A3, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects())


class SubStateNoOrdReacTestCase(ComplexReacTest5Subunits):
    """Test reaction declaration with SubUnitStates and no ordering."""
    def testFullyDeterminedFwdReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B, C = self.A, self.B, self.C
        with self.mdl:
            with self.vsys:
                with C[:, A2, A3, B1, B2]:
                    A1 >r[1]> A2
                r[1].K = 1
            expectedReacs = [
                ([C[A1, A2, A3, B1, B2]], [C[A2, A2, A3, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects())

            with self.vsys:
                with C[..., B1, B2]:
                    2*A1 + A3 >r[2]> 3*A2
                r[2].K = 1
            expectedReacs = [
                ([C[A1, A1, A3, B1, B2]], [C[A2, A2, A2, B1, B2]])
            ]
            self._checkReactions(expectedReacs, r[2]._getStepsObjects())

            with self.vsys:
                with self.assertRaises(Exception):
                    4*A1 >r[3]> A1 + 3*A2

    def testFullyDeterminedBidirReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3
        A, B, C = self.A, self.B, self.C
        with self.mdl:
            with self.vsys:
                with C[:, A2, A3, B1, B2]:
                    A1 <r[1]> A2
                r[1].K = 1, 1
            expectedReacs = [
                ([C[A1, A2, A3, B1, B2]], [C[A2, A2, A3, B1, B2]], 1),
                ([C[A1, A2, A3, B1, B2]], [C[A1, A1, A3, B1, B2]], 1),
                ([C[A2, A2, A3, B1, B2]], [C[A1, A2, A3, B1, B2]], 2),
                ([C[A3, A2, A3, B1, B2]], [C[A3, A1, A3, B1, B2]], 1),
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects(), checkRates=True)


    def testCompDepRate(self):
        r = self.r
        A1, A2, A3, B1, B2, B3, E1, E2 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3, self.E1, self.E2
        A, B, C, F = self.A, self.B, self.C, self.F
        S1, S2 = self.S1, self.S2
        with self.mdl:
            with self.vsys:
    
                with C[A3, ..., B1, B1] as C1, C[A3, ..., B2, B2] as C2:
                    A1[C1] + A1[C2] >r[1]> A2[C1] + A2[C2]
                    r[1].K = CompDepRate(lambda s: s.Count(A1), [C1])
            expectedReacs = [
                ([C[A3, A1, A1, B1, B1], C[A3, A1, A1, B2, B2]], [C[A3, A1, A2, B1, B1], C[A3, A1, A2, B2, B2]], 2*4),
                ([C[A3, A1, A1, B1, B1], C[A3, A1, A2, B2, B2]], [C[A3, A1, A2, B1, B1], C[A3, A2, A2, B2, B2]], 2*2),
                ([C[A3, A1, A1, B1, B1], C[A3, A1, A3, B2, B2]], [C[A3, A1, A2, B1, B1], C[A3, A2, A3, B2, B2]], 2*2),
                ([C[A3, A1, A2, B1, B1], C[A3, A1, A1, B2, B2]], [C[A3, A2, A2, B1, B1], C[A3, A2, A1, B2, B2]], 1*2),
                ([C[A3, A1, A2, B1, B1], C[A3, A1, A2, B2, B2]], [C[A3, A2, A2, B1, B1], C[A3, A2, A2, B2, B2]], 1*1),
                ([C[A3, A1, A2, B1, B1], C[A3, A1, A3, B2, B2]], [C[A3, A2, A2, B1, B1], C[A3, A2, A3, B2, B2]], 1*1),
                ([C[A3, A1, A3, B1, B1], C[A3, A1, A1, B2, B2]], [C[A3, A2, A3, B1, B1], C[A3, A2, A1, B2, B2]], 1*2),
                ([C[A3, A1, A3, B1, B1], C[A3, A1, A2, B2, B2]], [C[A3, A2, A3, B1, B1], C[A3, A2, A2, B2, B2]], 1*1),
                ([C[A3, A1, A3, B1, B1], C[A3, A1, A3, B2, B2]], [C[A3, A2, A3, B1, B1], C[A3, A2, A3, B2, B2]], 1*1),
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[A3, ..., B1, B1] as C1, C[A2, ..., B2, B2] as C2:
                    A1[C1] + A1[C2] >r[2]> A2[C1] + A2[C2]
                    r[2].K = CompDepRate(lambda s1, s2: s1.Count(A1) * s2.Count(A2), [C1, C2])
            expectedReacs = [
                ([C[A3, A1, A1, B1, B1], C[A2, A1, A1, B2, B2]], [C[A3, A1, A2, B1, B1], C[A2, A1, A2, B2, B2]], 2*1*4),
                ([C[A3, A1, A1, B1, B1], C[A2, A1, A2, B2, B2]], [C[A3, A1, A2, B1, B1], C[A2, A2, A2, B2, B2]], 2*2*2),
                ([C[A3, A1, A1, B1, B1], C[A2, A1, A3, B2, B2]], [C[A3, A1, A2, B1, B1], C[A2, A2, A3, B2, B2]], 2*1*2),
                ([C[A3, A1, A2, B1, B1], C[A2, A1, A1, B2, B2]], [C[A3, A2, A2, B1, B1], C[A2, A2, A1, B2, B2]], 1*1*2),
                ([C[A3, A1, A2, B1, B1], C[A2, A1, A2, B2, B2]], [C[A3, A2, A2, B1, B1], C[A2, A2, A2, B2, B2]], 1*2*1),
                ([C[A3, A1, A2, B1, B1], C[A2, A1, A3, B2, B2]], [C[A3, A2, A2, B1, B1], C[A2, A2, A3, B2, B2]], 1*1*1),
                ([C[A3, A1, A3, B1, B1], C[A2, A1, A1, B2, B2]], [C[A3, A2, A3, B1, B1], C[A2, A2, A1, B2, B2]], 1*1*2),
                ([C[A3, A1, A3, B1, B1], C[A2, A1, A2, B2, B2]], [C[A3, A2, A3, B1, B1], C[A2, A2, A2, B2, B2]], 1*2*1),
                ([C[A3, A1, A3, B1, B1], C[A2, A1, A3, B2, B2]], [C[A3, A2, A3, B1, B1], C[A2, A2, A3, B2, B2]], 1*1*1),
            ]
            self._checkReactions(expectedReacs, r[2]._getStepsObjects(), checkRates=True)


    def testPartiallyDeterminedFwdReac(self):
        r = self.r
        A1, A2, A3, B1, B2, B3, E1, E2 = self.A1, self.A2, self.A3, self.B1, self.B2, self.B3, self.E1, self.E2
        A, B, C, F = self.A, self.B, self.C, self.F
        S1, S2 = self.S1, self.S2
        with self.mdl:
            with self.vsys:
                with C[:, A2, :, B1, B2]:
                    A1 >r[1]> A2
                r[1].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B2]], [C[A2, A2, A1, B1, B2]], 2),
                ([C[A1, A2, A2, B1, B2]], [C[A2, A2, A2, B1, B2]], 1),
                ([C[A1, A2, A3, B1, B2]], [C[A2, A2, A3, B1, B2]], 1)
            ]
            self._checkReactions(expectedReacs, r[1]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[:, A2, A3, :, B2]:
                    A1 + ~B1 >r[2]> A2 + B1
                r[2].K = 1
            expectedReacs = [
                ([C[A1, A2, A3, B1, B2]], [C[A2, A2, A3, B1, B1]], 1),
                ([C[A1, A2, A3, B2, B2]], [C[A2, A2, A3, B1, B2]], 2),
                ([C[A1, A2, A3, B3, B2]], [C[A2, A2, A3, B1, B2]], 1),
                ([C[A1, A2, A3, B3, B2]], [C[A2, A2, A3, B3, B1]], 1),
            ]
            self._checkReactions(expectedReacs, r[2]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[:, A2, ..., B2]:
                    A1 + (B1|B3) >r[3]> A2 + B2
                r[3].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B2]], [C[A2, A2, A1, B2, B2]], 2),
                ([C[A1, A2, A2, B1, B2]], [C[A2, A2, A2, B2, B2]], 1),
                ([C[A1, A2, A3, B1, B2]], [C[A2, A2, A3, B2, B2]], 1),
                ([C[A1, A2, A1, B3, B2]], [C[A2, A2, A1, B2, B2]], 2),
                ([C[A1, A2, A2, B3, B2]], [C[A2, A2, A2, B2, B2]], 1),
                ([C[A1, A2, A3, B3, B2]], [C[A2, A2, A3, B2, B2]], 1)
            ]
            self._checkReactions(expectedReacs, r[3]._getStepsObjects(), checkRates=True)

            with self.vsys:
                C2 = C.get()
                with C[:, A2, A3, :, B2] as C3:
                    C2[A2, A2, ..., B2] + A1 >r[17]> C2[A1, A2, ..., B2] + A2
                r[17].K = 1
            expectedReacs = [
                ([C[A2, A2, A1, B1, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A1, B1, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A1, B1, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A1, B1, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A1, B1, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A1, B1, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A1, B2, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A1, B2, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A1, B2, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A1, B2, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A1, B2, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A1, B2, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A1, B3, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A1, B3, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A1, B3, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A1, B3, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A1, B3, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A1, B3, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A2, B1, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A2, B1, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A2, B1, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A2, B1, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A2, B1, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A2, B1, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A2, B2, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A2, B2, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A2, B2, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A2, B2, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A2, B2, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A2, B2, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A2, B3, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A2, B3, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A2, B3, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A2, B3, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A2, B3, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A2, B3, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A3, B1, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A3, B1, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A3, B1, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A3, B1, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A3, B1, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A3, B1, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A3, B2, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A3, B2, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A3, B2, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A3, B2, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A3, B2, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A3, B2, B2], C[A2, A2, A3, B3, B2]], 1),
                ([C[A2, A2, A3, B3, B2], C[A1, A2, A3, B1, B2]], [C[A1, A2, A3, B3, B2], C[A2, A2, A3, B1, B2]], 1),
                ([C[A2, A2, A3, B3, B2], C[A1, A2, A3, B2, B2]], [C[A1, A2, A3, B3, B2], C[A2, A2, A3, B2, B2]], 1),
                ([C[A2, A2, A3, B3, B2], C[A1, A2, A3, B3, B2]], [C[A1, A2, A3, B3, B2], C[A2, A2, A3, B3, B2]], 1),
            ]
            self._checkReactions(expectedReacs, r[17]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[A1, :, :, B1, B1]:
                    A1 >r[19]> A2
                r[19].K = 1
            expectedReacs = [
                ([C[A1, A1, A1, B1, B1]], [C[A1, A2, A1, B1, B1]], 3),
                ([C[A1, A1, A2, B1, B1]], [C[A1, A2, A2, B1, B1]], 2),
                ([C[A1, A1, A3, B1, B1]], [C[A1, A2, A3, B1, B1]], 2),
                ([C[A1, A2, A2, B1, B1]], [C[A2, A2, A2, B1, B1]], 1),
                ([C[A1, A2, A3, B1, B1]], [C[A2, A2, A3, B1, B1]], 1),
                ([C[A1, A3, A3, B1, B1]], [C[A2, A3, A3, B1, B1]], 1),
            ]
            self._checkReactions(expectedReacs, r[19]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[...]:
                    2*A1 + (B1|B3) >r[20]> 2*A2 + B2
                r[20].K = 1
            expectedReacs = [
                ([C[A1, A1, A1, B1, B1]], [C[A2, A2, A1, B2, B1]], 6),
                ([C[A1, A1, A2, B1, B1]], [C[A2, A2, A2, B2, B1]], 2),
                ([C[A1, A1, A3, B1, B1]], [C[A2, A2, A3, B2, B1]], 2),
                ([C[A1, A1, A1, B1, B2]], [C[A2, A2, A1, B2, B2]], 3),
                ([C[A1, A1, A2, B1, B2]], [C[A2, A2, A2, B2, B2]], 1),
                ([C[A1, A1, A3, B1, B2]], [C[A2, A2, A3, B2, B2]], 1),
                ([C[A1, A1, A1, B1, B3]], [C[A2, A2, A1, B2, B1]], 3),
                ([C[A1, A1, A1, B1, B3]], [C[A2, A2, A1, B2, B3]], 3),
                ([C[A1, A1, A2, B1, B3]], [C[A2, A2, A2, B2, B1]], 1),
                ([C[A1, A1, A2, B1, B3]], [C[A2, A2, A2, B2, B3]], 1),
                ([C[A1, A1, A3, B1, B3]], [C[A2, A2, A3, B2, B1]], 1),
                ([C[A1, A1, A3, B1, B3]], [C[A2, A2, A3, B2, B3]], 1),
                ([C[A1, A1, A1, B2, B3]], [C[A2, A2, A1, B2, B2]], 3),
                ([C[A1, A1, A2, B2, B3]], [C[A2, A2, A2, B2, B2]], 1),
                ([C[A1, A1, A3, B2, B3]], [C[A2, A2, A3, B2, B2]], 1),
                ([C[A1, A1, A1, B3, B3]], [C[A2, A2, A1, B2, B3]], 6),
                ([C[A1, A1, A2, B3, B3]], [C[A2, A2, A2, B2, B3]], 2),
                ([C[A1, A1, A3, B3, B3]], [C[A2, A2, A3, B2, B3]], 2),
            ]
            self._checkReactions(expectedReacs, r[20]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[..., B1, B1]:
                    ~A1 >r[21]> A1
                r[21].K = 1
            expectedReacs = [
                ([C[A2, A1, A1, B1, B1]], [C[A1, A1, A1, B1, B1]], 1),
                ([C[A3, A1, A1, B1, B1]], [C[A1, A1, A1, B1, B1]], 1),
                ([C[A2, A2, A1, B1, B1]], [C[A2, A1, A1, B1, B1]], 2),
                ([C[A2, A3, A1, B1, B1]], [C[A2, A1, A1, B1, B1]], 1),
                ([C[A2, A3, A1, B1, B1]], [C[A3, A1, A1, B1, B1]], 1),
                ([C[A3, A3, A1, B1, B1]], [C[A3, A1, A1, B1, B1]], 2),
                ([C[A2, A2, A2, B1, B1]], [C[A2, A2, A1, B1, B1]], 3),
                ([C[A2, A2, A3, B1, B1]], [C[A2, A3, A1, B1, B1]], 2),
                ([C[A2, A2, A3, B1, B1]], [C[A2, A2, A1, B1, B1]], 1),
                ([C[A2, A3, A3, B1, B1]], [C[A2, A3, A1, B1, B1]], 2),
                ([C[A2, A3, A3, B1, B1]], [C[A3, A3, A1, B1, B1]], 1),
                ([C[A3, A3, A3, B1, B1]], [C[A3, A3, A1, B1, B1]], 3),
            ]
            self._checkReactions(expectedReacs, r[21]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[A1,...,B1,B1]:
                    A1 + (A1|A2) >r[22]> A3 + (A1|A2)
                r[22].K = 1
            expectedReacs = [
                ([C[A1, A1, A1, B1, B1]], [C[A3, A1, A1, B1, B1]], 3),
                ([C[A1, A1, A2, B1, B1]], [C[A3, A1, A2, B1, B1]], 3),
                ([C[A1, A1, A3, B1, B1]], [C[A3, A1, A3, B1, B1]], 1),
                ([C[A1, A2, A2, B1, B1]], [C[A3, A2, A2, B1, B1]], 2),
                ([C[A1, A2, A3, B1, B1]], [C[A3, A2, A3, B1, B1]], 1),
            ]
            self._checkReactions(expectedReacs, r[22]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[...,B1,B1]:
                    A1 + S1 >r[23]> A1 + S2
                r[23].K = 1
            expectedReacs = [
                ([C[A1, A1, A1, B1, B1], S1], [C[A1, A1, A1, B1, B1], S2], 3),
                ([C[A1, A1, A2, B1, B1], S1], [C[A1, A1, A2, B1, B1], S2], 2),
                ([C[A1, A1, A3, B1, B1], S1], [C[A1, A1, A3, B1, B1], S2], 2),
                ([C[A1, A2, A2, B1, B1], S1], [C[A1, A2, A2, B1, B1], S2], 1),
                ([C[A1, A2, A3, B1, B1], S1], [C[A1, A2, A3, B1, B1], S2], 1),
                ([C[A1, A3, A3, B1, B1], S1], [C[A1, A3, A3, B1, B1], S2], 1),
            ]
            self._checkReactions(expectedReacs, r[23]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[:, :, :, B1, B1]:
                    A1['a'] + A2['b'] >r[24]> A2['a'] + A3['b']
                r[24].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B1]], [C[A2, A3, A1, B1, B1]], 2),
                ([C[A1, A2, A2, B1, B1]], [C[A2, A3, A2, B1, B1]], 2),
                ([C[A1, A2, A3, B1, B1]], [C[A2, A3, A3, B1, B1]], 1),
            ]
            self._checkReactions(expectedReacs, r[24]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[:, :, :, B1, B1]:
                    A1['a'] + A2['b'] >r[57]> A3['a'] + A3['b']
                    A3['a'] + A3['b'] >r[58]> A1['a'] + A2['b']
                r[57].K = 1
                r[58].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B1]], [C[A3, A3, A1, B1, B1]], 2),
                ([C[A1, A2, A2, B1, B1]], [C[A3, A3, A2, B1, B1]], 2),
                ([C[A1, A2, A3, B1, B1]], [C[A3, A3, A3, B1, B1]], 1),
            ]
            self._checkReactions(expectedReacs, r[57]._getStepsObjects(), checkRates=True)
            expectedReacs = [
                ([C[A3, A3, A1, B1, B1]], [C[A1, A2, A1, B1, B1]], 1),
                ([C[A3, A3, A2, B1, B1]], [C[A1, A2, A2, B1, B1]], 1),
                ([C[A3, A3, A3, B1, B1]], [C[A1, A2, A3, B1, B1]], 3),
            ]
            self._checkReactions(expectedReacs, r[58]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[...]:
                    A1 + B1 >r[25]> A2 + B2
                r[25].K = 1
            expectedReacs = [
                ([C[A1, A1, A1, B1, B1]], [C[A2, A1, A1, B2, B1]], 6),
                ([C[A1, A1, A1, B1, B2]], [C[A2, A1, A1, B2, B2]], 3),
                ([C[A1, A1, A1, B1, B3]], [C[A2, A1, A1, B2, B3]], 3),
                ([C[A1, A1, A2, B1, B1]], [C[A2, A1, A2, B2, B1]], 4),
                ([C[A1, A1, A2, B1, B2]], [C[A2, A1, A2, B2, B2]], 2),
                ([C[A1, A1, A2, B1, B3]], [C[A2, A1, A2, B2, B3]], 2),
                ([C[A1, A1, A3, B1, B1]], [C[A2, A1, A3, B2, B1]], 4),
                ([C[A1, A1, A3, B1, B2]], [C[A2, A1, A3, B2, B2]], 2),
                ([C[A1, A1, A3, B1, B3]], [C[A2, A1, A3, B2, B3]], 2),
                ([C[A1, A2, A2, B1, B1]], [C[A2, A2, A2, B2, B1]], 2),
                ([C[A1, A2, A2, B1, B2]], [C[A2, A2, A2, B2, B2]], 1),
                ([C[A1, A2, A2, B1, B3]], [C[A2, A2, A2, B2, B3]], 1),
                ([C[A1, A2, A3, B1, B1]], [C[A2, A2, A3, B2, B1]], 2),
                ([C[A1, A2, A3, B1, B2]], [C[A2, A2, A3, B2, B2]], 1),
                ([C[A1, A2, A3, B1, B3]], [C[A2, A2, A3, B2, B3]], 1),
                ([C[A1, A3, A3, B1, B1]], [C[A2, A3, A3, B2, B1]], 2),
                ([C[A1, A3, A3, B1, B2]], [C[A2, A3, A3, B2, B2]], 1),
                ([C[A1, A3, A3, B1, B3]], [C[A2, A3, A3, B2, B3]], 1),
            ]
            self._checkReactions(expectedReacs, r[25]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[A3, ..., B1, B1] as C1, C[A3, ..., B1, B1] as C2:
                    A1[C1] + A1[C2] >r[26]> A2[C1] + A2[C2]
                r[26].K = 1
            expectedReacs = [
                ([C[A3, A1, A1, B1, B1], C[A3, A1, A1, B1, B1]], [C[A3, A1, A2, B1, B1], C[A3, A1, A2, B1, B1]], 2),
                ([C[A3, A1, A1, B1, B1], C[A3, A1, A2, B1, B1]], [C[A3, A1, A2, B1, B1], C[A3, A2, A2, B1, B1]], 2),
                ([C[A3, A1, A1, B1, B1], C[A3, A1, A3, B1, B1]], [C[A3, A1, A2, B1, B1], C[A3, A2, A3, B1, B1]], 2),
                ([C[A3, A1, A2, B1, B1], C[A3, A1, A2, B1, B1]], [C[A3, A2, A2, B1, B1], C[A3, A2, A2, B1, B1]], 0.5),
                ([C[A3, A1, A2, B1, B1], C[A3, A1, A3, B1, B1]], [C[A3, A2, A2, B1, B1], C[A3, A2, A3, B1, B1]], 1),
                ([C[A3, A1, A3, B1, B1], C[A3, A1, A3, B1, B1]], [C[A3, A2, A3, B1, B1], C[A3, A2, A3, B1, B1]], 0.5),
            ]
            self._checkReactions(expectedReacs, r[26]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with C[:, ~A1, ~A1, B1, B1]:
                    A2 >r[70]> A3
                r[70].K = 1
            expectedReacs = [
                ([C[A1, A2, A3, B1, B1]], [C[A1, A3, A3, B1, B1]], 1),
                ([C[A1, A2, A2, B1, B1]], [C[A1, A2, A3, B1, B1]], 2),
                ([C[A2, A3, A3, B1, B1]], [C[A3, A3, A3, B1, B1]], 1),
                ([C[A2, A2, A3, B1, B1]], [C[A2, A3, A3, B1, B1]], 2),
                ([C[A2, A2, A2, B1, B1]], [C[A2, A2, A3, B1, B1]], 3),
            ]
            self._checkReactions(expectedReacs, r[70]._getStepsObjects(), checkRates=True)

            with self.vsys:
                with self.assertRaises(Exception):
                    with C[A1, A2, ..., B2]:
                        A1 + B1|B3 >r[4]> A1 + B2

                with self.assertRaises(Exception):
                    with C[A1, A2, ..., B2]:
                        B1|B3 + A1 >r[40]> A1 + B2

                with self.assertRaises(Exception):
                    with C[...]:
                        2*A1 >r[41]> A2

                with self.assertRaises(Exception):
                    with C[...]:
                        A2 >r[42]> 2*A1

                with self.assertRaises(Exception):
                    with C[A1, A2, ..., B2]:
                        B1 >r[5]> A1

                with self.assertRaises(Exception):
                    with C[A1, A2, ..., B2]:
                        B1 >r[6]> None

                with self.assertRaises(Exception):
                    with C[A1, A2, ..., B2]:
                        B1|B3 >r[7]> ~B1

                with self.assertRaises(Exception):
                    with C[:, A2, :, B1, B2] as C2, C[:, A2, :, B1, B2] as C3:
                        A1[C2] >r[8]> A2[C3]

                with self.assertRaises(Exception):
                    A1 >r[9]> A2

                with self.assertRaises(Exception):
                    with C[A1, A2, ..., B2]:
                        E1 >r[10]> E2

                with self.assertRaises(Exception):
                    with C[:, A2, :, B1, B2] as C2, C[:, A2, :, B1, B2] as C3:
                        A1 >r[11]> A2

                with self.assertRaises(Exception):
                    C2[:, A2, ..., B2] + A1 >r[15]> C2[A1, A2, ..., B2]

                C2 = C.get()
                with self.assertRaises(Exception):
                    C2[:, A2, ..., B2] + A1[C2] >r[16]> C2[A1, A2, ..., B2]

                C2 = C.get()
                with self.assertRaises(Exception):
                    C2[A1, A2, A3, ...] + C2[A3, A3, A3, ...] >r[18]> C2[A2, A2, A2, ...] + C2[A1, A1, A1, ...]

                with self.assertRaises(Exception):
                    with C[~A1, ~A1, ~A1, B1, B1]:
                        A1 >r[50]> A2
                        r[50].K = 1

            with self.ssys:
                with C[:, A2, :, B1, B2]:
                    A1.s >r[13]> A2.s
                r[13].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B2].s], [C[A2, A2, A1, B1, B2].s], 2),
                ([C[A1, A2, A2, B1, B2].s], [C[A2, A2, A2, B1, B2].s], 1),
                ([C[A1, A2, A3, B1, B2].s], [C[A2, A2, A3, B1, B2].s], 1)
            ]
            self._checkSurfReactions(expectedReacs, r[13]._getStepsObjects(), checkRates=True)

            with self.ssys:
                with C[:, A2, :, B1, B2]:
                    A1.i >r[14]> A2.o
                r[14].K = 1
            expectedReacs = [
                ([C[A1, A2, A1, B1, B2].i], [C[A2, A2, A1, B1, B2].o], 2),
                ([C[A1, A2, A2, B1, B2].i], [C[A2, A2, A2, B1, B2].o], 1),
                ([C[A1, A2, A3, B1, B2].i], [C[A2, A2, A3, B1, B2].o], 1)
            ]
            self._checkSurfReactions(expectedReacs, r[14]._getStepsObjects(), checkRates=True)

            with self.ssys:
                with self.assertRaises(Exception):
                    with C[:, A2, :, B1, B2]:
                        A1.s + A1.o >r[12]> A2.s + A1.s



def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ComplexCreation, "test"))
    all_tests.append(unittest.makeSuite(ComplexCreationStatesAsSpecies, "test"))
    all_tests.append(unittest.makeSuite(CompStateGeneral, "test"))
    all_tests.append(unittest.makeSuite(CompStateNoOrdering, "test"))
    all_tests.append(unittest.makeSuite(CompStateStrongOrdering, "test"))
    all_tests.append(unittest.makeSuite(CompStateRotationalSymmetryOrdering, "test"))
    all_tests.append(unittest.makeSuite(CompSelNoOrdReacTestCase, "test"))
    all_tests.append(unittest.makeSuite(SubStateNoOrdReacTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
