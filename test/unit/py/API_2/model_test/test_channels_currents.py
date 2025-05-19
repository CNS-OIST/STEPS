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

""" Unit tests for declaration and simulation of channels and currents."""

import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.utils import *

import steps.API_1.model as smodel

class ChannelCurrentDeclarations(unittest.TestCase):
    """Test model checking."""
    def setUp(self):
        pass

    def testChannelDeclaration(self):
        mdl = Model()

        with mdl:
            sus1, sus2, sus3, sus4 = SubUnitState.Create()
            suA, suB = SubUnit.Create([sus1, sus2], [sus3, sus4])

            with self.assertRaises(NotImplementedError):
                CC = Channel.Create([suA, suB], statesAsSpecies=False)

            Chan = Channel.Create([suA, suB])

            # Load from STEPS object
            smdl = mdl._getStepsObjects()[0]
            chan2 = smodel.Chan('chan2', smdl)
            chanstate1 = smodel.ChanState('chanstate1', smdl, chan2)
            chanstate2 = smodel.ChanState('chanstate2', smdl, chan2)
            chanstate3 = smodel.ChanState('chanstate3', smdl, chan2)

            chan2b = Channel._FromStepsObject(chan2, mdl)
            self.assertEqual(chan2b.name, 'chan2')
            states = [state._getStepsObjects()[0] for state in chan2b[...]]
            stepsStates = [chanstate1, chanstate2, chanstate3]
            self.assertCountEqual(states, stepsStates)

    def testVDepReacDeclaration(self):
        mdl = Model()
        r = ReactionManager()

        with mdl:
            SA, SB = Species.Create()

            vsys = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

            with ssys:
                SA.s >r[1]> SB.s
                r[1].K = VDepRate(lambda V: V * 3)

                self.assertEqual(r[1].K(10), 30)

                SA.s <r[2]> SB.s
                r[2].K = VDepRate(lambda V: V * 2), VDepRate(lambda V: V * 4)

                va, vb = r[2].K
                self.assertEqual(va(10), 20)
                self.assertEqual(vb(10), 40)

            with self.assertRaises(NotImplementedError):
                with vsys:
                    SA >r[3]> SB
                    r[3].K = VDepRate(lambda V: V * 3)

    def testOhmicCurrentDeclaration(self):
        mdl = Model()

        with mdl:
            sus1, sus2, sus3, sus4, sus5, sus6, sus7, sus8 = SubUnitState.Create()
            suA, suB, suC, suD = SubUnit.Create([sus1, sus2], [sus3, sus4], [sus5, sus6], [sus7, sus8])

            Chan1, Chan2 = Channel.Create([suA, suB], [suC, suD])

            ssys = SurfaceSystem.Create()

            with ssys:
                Curr1 = OhmicCurr.Create(Chan1[sus1, sus3], 1e-10, -70e-3)
            self.assertEqual(Curr1[sus1, sus3]._getStepsObjects()[0].g, 1e-10)
            self.assertEqual(Curr1[sus1, sus3]._getStepsObjects()[0].erev, -70e-3)
            with self.assertRaises(Exception):
                Curr1[sus1, sus4]

            with self.assertRaises(TypeError):
                with ssys:
                    OhmicCurr(Chan1[sus1, sus3], 'value42', -75e-3)

            with self.assertRaises(TypeError):
                with ssys:
                    OhmicCurr(Chan1[sus1, sus3], 1e-10, 'value')

            # Default values
            with ssys:
                OhmicCurr(Chan1[sus1, sus3], 1e-10)
            with ssys:
                OhmicCurr(Chan1[sus1, sus3])

            cdc1 = CompDepCond(lambda s: s.Count(sus1) + 10*s.Count(sus3), [Chan1])
            with self.assertRaises(TypeError):
                CompDepCond(42, [Chan1])
            with self.assertRaises(TypeError):
                CompDepCond(lambda s: 1, 42)

            with self.assertRaises(TypeError):
                with ssys:
                    cdc2 = CompDepCond(lambda s: 1, [42])
                    OhmicCurr(Chan1[sus1, sus3], cdc2, -70e-3)

            with ssys:
                Curr2 = OhmicCurr(Chan1, cdc1, -70e-3)
            self.assertEqual(Curr2[sus1, sus3]._getStepsObjects()[0].g, 11)
            self.assertEqual(Curr2[sus1, sus4]._getStepsObjects()[0].g, 1)
            self.assertEqual(Curr2[sus2, sus3]._getStepsObjects()[0].g, 10)
            self.assertEqual(Curr2[sus2, sus4]._getStepsObjects()[0].g, 0)
            self.assertEqual(Curr2[sus1, sus3]._getStepsObjects()[0].erev, -70e-3)
            self.assertEqual(Curr2[sus1, sus4]._getStepsObjects()[0].erev, -70e-3)
            self.assertEqual(Curr2[sus2, sus3]._getStepsObjects()[0].erev, -70e-3)
            self.assertEqual(Curr2[sus2, sus4]._getStepsObjects()[0].erev, -70e-3)
            self.assertEqual(set(curr._getStepsObjects()[0].g for curr in Curr2[sus1,:]), set([11,1]))
            self.assertEqual(set(curr._getStepsObjects()[0].g for curr in Curr2), set([11,1,10,0]))

            with ssys:
                Curr3 = OhmicCurr(Chan1[sus1,:], cdc1, -70e-3)
            self.assertEqual(Curr2[sus1, sus3]._getStepsObjects()[0].g, 11)
            self.assertEqual(Curr2[sus1, sus4]._getStepsObjects()[0].g, 1)
            self.assertEqual(Curr2[sus1, sus3]._getStepsObjects()[0].erev, -70e-3)
            self.assertEqual(Curr2[sus1, sus4]._getStepsObjects()[0].erev, -70e-3)
            with self.assertRaises(Exception):
                Curr3[sus2, sus3]
            with self.assertRaises(Exception):
                Curr3[sus2, sus4]

            with self.assertRaises(Exception):
                with ssys:
                    OhmicCurr(Chan2, cdc1, -70e-3)

            with self.assertRaises(TypeError):
                with ssys:
                    OhmicCurr(42, 1e-10, -70e-3)

    def testGHCKurrentDeclaration(self):
        mdl = Model()

        with mdl:
            Ca, noVal = Species.Create()
        
            with self.assertRaises(TypeError):
                Ca.valence = 'test'
            Ca.valence = 2

            self.assertEqual(Ca.valence, 2)

            Na = Species.Create(valence=1)
            self.assertEqual(Na.valence, 1)

            sus1, sus2, sus3, sus4, sus5, sus6, sus7, sus8 = SubUnitState.Create()
            suA, suB, suC, suD = SubUnit.Create([sus1, sus2], [sus3, sus4], [sus5, sus6], [sus7, sus8])

            Chan1, Chan2 = Channel.Create([suA, suB], [suC, suD])

            ssys = SurfaceSystem.Create()

            with self.assertRaises(Exception):
                with ssys:
                    GHKCurr(Chan1[sus1, sus3], noVal, 1e-20)

            with ssys:
                Curr1 = GHKCurr.Create(Chan1[sus1, sus3], Ca, 1e-20)
            with self.assertRaises(Exception):
                Curr1[sus1, sus4]

            with self.assertRaises(TypeError):
                with ssys:
                    GHKCurr(Chan1[sus1, sus3], 'value', 1e-20)
            with self.assertRaises(TypeError):
                with ssys:
                    GHKCurr(Chan1[sus1, sus3], Ca, 'value')
            with self.assertRaises(Exception):
                with ssys:
                    GHKCurr(Chan1[sus1, sus3], Ca)
            with self.assertRaises(Exception):
                with ssys:
                    GHKCurr(Chan1[sus1, sus3])

            cdc1 = CompDepP(lambda s: 1 + s.Count(sus1) + 10*s.Count(sus3), [Chan1])
            with self.assertRaises(TypeError):
                CompDepP(42, [Chan1])
            with self.assertRaises(TypeError):
                CompDepP(lambda s: 1, 42)

            with self.assertRaises(TypeError):
                cdc2 = CompDepP(lambda s: 1, [42])
                with ssys:
                    GHKCurr(Chan1[sus1, sus3], Ca, cdc2)

            with ssys:
                Curr2 = GHKCurr(Chan1, Ca, cdc1)

                Curr3 = GHKCurr(Chan1[sus1,:], Ca, cdc1)

            with self.assertRaises(Exception):
                Curr3[sus2, sus3]
            with self.assertRaises(Exception):
                Curr3[sus2, sus4]

            with self.assertRaises(Exception):
                with ssys:
                    GHKCurr(Chan2, Ca, cdc1)

            with self.assertRaises(TypeError):
                with ssys:
                    GHKCurr(42, Ca, 1e-20)

            with ssys:
                infos = GHKCurr.PInfo(
                    g = 20e-12, V = -22e-3, T = 293.15, oconc = 4e-3, iconc = 155e-3
                )
                # Cannot access the value since infos is not associated with a current
                with self.assertRaises(Exception):
                    infos.value

                Curr4 = GHKCurr(
                    Chan1[sus1, sus3],
                    Ca,
                    infos,
                )
                self.assertEqual(Curr4.P.value, 3.3006027395421925e-20)

                # Reuse the same infos to create another current
                Curr5 = GHKCurr(
                    Chan1[sus2, sus3],
                    Na,
                    infos,
                )
                self.assertEqual(Curr5.P.value, 9.009267378478541e-20)

            self.assertEqual(Curr4.P.value, Curr4._getStepsObjects()[0].getP())
            self.assertEqual(Curr5.P.value, Curr5._getStepsObjects()[0].getP())

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(ChannelCurrentDeclarations))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


