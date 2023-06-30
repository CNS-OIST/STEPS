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

""" Unit tests for vesicle objects declaration."""

import unittest

from steps import interface

from steps.model import *


class VesicleObjectsDeclaration(unittest.TestCase):
    """Test vesicle objects declaration."""
    def setUp(self):
        self.mdl = Model()

    def testVesicleCreation(self):
        with self.mdl:
            vssys1 = VesicleSurfaceSystem.Create()

            ves1 = Vesicle.Create(1, 2)
            self.assertEqual(ves1.Diameter, 1)
            self.assertEqual(ves1.Dcst, 2)
            self.assertEqual(set(ves1.sysNames), set())
            ves1.addSystem(vssys1)
            self.assertEqual(set(ves1._getStepsObjects()[0].getVesSurfsys()), set([vssys1.name]))

            ves1.Dcst = 3
            self.assertEqual(ves1.Dcst, 3)

            with self.assertRaises(TypeError):
                ves1.addSystem(42)

            with self.assertRaises(TypeError):
                ves1.Dcst = '42'

            ves2 = Vesicle.Create(2, 3, vssys1)
            self.assertEqual(ves2.Diameter, 2)
            self.assertEqual(ves2.Dcst, 3)

            stepsVes = ves1._getStepsObjects()[0].__class__('ves1b', self.mdl._getStepsObjects()[0], ves1.Diameter, ves1.Dcst)
            stepsVes.addVesSurfsys('vssys1')

            ves1b = Vesicle._FromStepsObject(stepsVes, self.mdl)
            self.assertEqual(ves1b.Diameter, ves1.Diameter)
            self.assertEqual(ves1b.Dcst, ves1.Dcst)
            self.assertEqual(set(ves1b.sysNames), set(ves1.sysNames))
            self.assertEqual(set(ves1b._getStepsObjects()[0].getVesSurfsys()), set([vssys1.name]))

        with self.assertRaises(Exception):
            ves = Vesicle.Create()


    def testRaftCreation(self):
        with self.mdl:
            rssys1 = RaftSurfaceSystem.Create()

            raft1 = Raft.Create(4, 5)
            self.assertEqual(raft1.Diameter, 4)
            self.assertEqual(raft1.Dcst, 5)
            raft1.addSystem(rssys1)
            self.assertEqual(set(raft1._getStepsObjects()[0].getRaftsys()), set([rssys1.name]))

            raft1.Dcst = 3
            self.assertEqual(raft1.Dcst, 3)

            with self.assertRaises(TypeError):
                raft1.addSystem(42)

            with self.assertRaises(TypeError):
                raft1.Dcst = '42'

            raft2 = Raft.Create(6, 7, rssys1)
            self.assertEqual(raft2.Diameter, 6)
            self.assertEqual(raft2.Dcst, 7)

            stepsRaft = raft1._getStepsObjects()[0].__class__('raft1b', self.mdl._getStepsObjects()[0], raft1.Diameter, raft1.Dcst)
            stepsRaft.addRaftsys('rssys1')

            raft1b = Raft._FromStepsObject(stepsRaft, self.mdl)
            self.assertEqual(raft1b.Diameter, raft1.Diameter)
            self.assertEqual(raft1b.Dcst, raft1.Dcst)
            self.assertEqual(set(raft1b.sysNames), set(raft1.sysNames))
            self.assertEqual(set(raft1b._getStepsObjects()[0].getRaftsys()), set([rssys1.name]))

        with self.assertRaises(Exception):
            raft = Raft.Create()

    def testEndocytosisCreation(self):
        with self.mdl:
            S1, S2 = Species.Create()
            ves1 = Vesicle.Create(1)
            ssys = SurfaceSystem.Create()

            with ssys:
                endo1 = Endocytosis.Create(ves1, 1, S1 + 2*S2, True)

                self.assertEqual(endo1.Vesicle.name, ves1.name)
                self.assertEqual(endo1.K, 1)
                self.assertEqual(sorted([spec.getID() for spec in endo1.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertTrue(endo1._getStepsObjects()[0].getInner())

                endo1.K = 2
                self.assertEqual(endo1.K, 2)

                with self.assertRaises(Exception):
                    endo1.Dependencies = S1 + S2

                with self.assertRaises(TypeError):
                    endo2 = Endocytosis.Create(42, 1, S1 + 2*S2)
                with self.assertRaises(TypeError):
                    endo3 = Endocytosis.Create(ves1, 'rate', S1 + 2*S2)
                with self.assertRaises(TypeError):
                    endo4 = Endocytosis.Create(ves1, 1, 'deps')
                with self.assertRaises(Exception):
                    endo5 = Endocytosis.Create(ves1, 1, S1.o + 2* S2.i, True)

                stepsEndo = endo1._getStepsObjects()[0].__class__('endo1b', ssys._getStepsObjects()[0], ves1._getStepsObjects()[0], None, endo1.Dependencies._GetStepsElems(), endo1.K)

                endo1b = Endocytosis._FromStepsObject(stepsEndo, self.mdl)
                self.assertEqual(endo1b.Vesicle.name, endo1.Vesicle.name)
                self.assertEqual(endo1b.K, endo1.K)
                self.assertEqual(sorted([spec.getID() for spec in endo1b.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertTrue(endo1b._getStepsObjects()[0].getInner())

            with self.assertRaises(Exception):
                endo = Endocytosis.Create(ves1, 1, S1 + 2*S2, True)

    def testExocytosisCreation(self):
        with self.mdl:
            S1, S2 = Species.Create()
            vssys = VesicleSurfaceSystem.Create()
            raft1 = Raft.Create(1)

            with vssys:
                exo1 = Exocytosis.Create(1, S1 + 2*S2, raft1)

                self.assertEqual(exo1.K, 1)
                self.assertEqual(sorted([spec.getID() for spec in exo1.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertEqual(exo1.Raft.name, raft1.name)

                exo1.K = 2
                self.assertEqual(exo1.K, 2)

                with self.assertRaises(Exception):
                    exo1.Dependencies = S1 + S2

                with self.assertRaises(TypeError):
                    exo2 = Exocytosis.Create('rate', S1 + 2*S2)
                with self.assertRaises(TypeError):
                    exo3 = Exocytosis.Create(3, 'deps')
                with self.assertRaises(TypeError):
                    exo4 = Exocytosis.Create(4, S1 + 2*S2, 'raft')
                with self.assertRaises(Exception):
                    exo5 = Exocytosis.Create(1, S1.o + 2* S2.i, raft1)

                stepsExo = exo1._getStepsObjects()[0].__class__('exo1b', vssys._getStepsObjects()[0], exo1.Dependencies._GetStepsElems(), raft1._getStepsObjects()[0], exo1.K)

                exo1b = Exocytosis._FromStepsObject(stepsExo, self.mdl)
                self.assertEqual(exo1b.Raft.name, exo1.Raft.name)
                self.assertEqual(exo1b.K, exo1.K)
                self.assertEqual(sorted([spec.getID() for spec in exo1b.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])

                # No raft
                exo6 = Exocytosis.Create(1, S1 + 2*S2)

                self.assertEqual(exo6.K, 1)
                self.assertEqual(sorted([spec.getID() for spec in exo6.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertEqual(exo6.Raft, None)

                exo6.K = 2
                self.assertEqual(exo6.K, 2)

                with self.assertRaises(Exception):
                    exo7.Dependencies = S2 + S2

                stepsExo = exo6._getStepsObjects()[0].__class__('exo6b', vssys._getStepsObjects()[0], exo6.Dependencies._GetStepsElems(), None, exo6.K)

                exo6b = Exocytosis._FromStepsObject(stepsExo, self.mdl)
                self.assertEqual(exo6b.Raft, None)
                self.assertEqual(exo6b.K, exo6.K)
                self.assertEqual(sorted([spec.getID() for spec in exo6b.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])

            with self.assertRaises(Exception):
                exo = Exocytosis.Create(1, S1 + 2*S2, raft1)

    def testRaftEndocytosisCreation(self):
        with self.mdl:
            S1, S2 = Species.Create()
            ves1 = Vesicle.Create(1)
            rssys = RaftSurfaceSystem.Create()

            with rssys:
                endo1 = RaftEndocytosis.Create(ves1, 1, S1 + 2*S2, True)

                self.assertEqual(endo1.Vesicle.name, ves1.name)
                self.assertEqual(endo1.K, 1)
                self.assertEqual(sorted([spec.getID() for spec in endo1.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertTrue(endo1._getStepsObjects()[0].getInner())

                endo1.K = 2
                self.assertEqual(endo1.K, 2)

                with self.assertRaises(Exception):
                    endo1.Dependencies = S1 + S2

                with self.assertRaises(TypeError):
                    endo2 = RaftEndocytosis.Create(42, 1, S1 + 2*S2)
                with self.assertRaises(TypeError):
                    endo3 = RaftEndocytosis.Create(ves1, 'rate', S1 + 2*S2)
                with self.assertRaises(TypeError):
                    endo4 = RaftEndocytosis.Create(ves1, 1, 'deps')
                with self.assertRaises(Exception):
                    endo5 = RaftEndocytosis.Create(ves1, 1, S1.o + 2* S2.i, True)

                stepsRaftEndo = endo1._getStepsObjects()[0].__class__('endo1b', rssys._getStepsObjects()[0], ves1._getStepsObjects()[0], None, endo1.Dependencies._GetStepsElems(), endo1.K)

                endo1b = RaftEndocytosis._FromStepsObject(stepsRaftEndo, self.mdl)
                self.assertEqual(endo1b.Vesicle.name, endo1.Vesicle.name)
                self.assertEqual(endo1b.K, endo1.K)
                self.assertEqual(sorted([spec.getID() for spec in endo1b.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertTrue(endo1b._getStepsObjects()[0].getInner())

            with self.assertRaises(Exception):
                endo = RaftEndocytosis.Create(ves1, 1, S1 + 2*S2, True)

    def testRaftGenCreation(self):
        with self.mdl:
            S1, S2 = Species.Create()
            vssys = SurfaceSystem.Create()
            raft1 = Raft.Create(1)

            with vssys:
                rgen1 = RaftGen.Create(raft1, 1, S1 + 2*S2)

                self.assertEqual(rgen1.K, 1)
                self.assertEqual(sorted([spec.getID() for spec in rgen1.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])
                self.assertEqual(rgen1.Raft.name, raft1.name)

                rgen1.K = 2
                self.assertEqual(rgen1.K, 2)

                with self.assertRaises(Exception):
                    rgen1.Dependencies = S1 + S2

                with self.assertRaises(TypeError):
                    rgen2 = RaftGen.Create('raft', 1, S1 + 2*S2)
                with self.assertRaises(TypeError):
                    rgen3 = RaftGen.Create(raft1, 'rate', S1 + 2*S2)
                with self.assertRaises(TypeError):
                    rgen4 = RaftGen.Create(raft1, 4, 'deps')
                with self.assertRaises(Exception):
                    rgen5 = RaftGen.Create(raft1, 1, S1.o + 2* S2.i)

                stepsExo = rgen1._getStepsObjects()[0].__class__('rgen1b', vssys._getStepsObjects()[0], rgen1.Dependencies._GetStepsElems(), raft1._getStepsObjects()[0], rgen1.K)

                rgen1b = RaftGen._FromStepsObject(stepsExo, self.mdl)
                self.assertEqual(rgen1b.Raft.name, rgen1.Raft.name)
                self.assertEqual(rgen1b.K, rgen1.K)
                self.assertEqual(sorted([spec.getID() for spec in rgen1b.Dependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])

            with self.assertRaises(Exception):
                rgen = RaftGen.Create(raft1, 1, S1 + 2*S2)

    def testRaftDisCreation(self):
        with self.mdl:
            S1, S2 = Species.Create()
            rssys = RaftSurfaceSystem.Create()
            raft1 = Raft.Create(1)

            with rssys:
                rdis1 = RaftDis.Create(1, S1 + 2*S2)

                self.assertEqual(rdis1.K, 1)
                self.assertEqual(sorted([spec.getID() for spec in rdis1.AntiDependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])

                rdis1.K = 2
                self.assertEqual(rdis1.K, 2)

                with self.assertRaises(Exception):
                    rdis1.AntiDependencies = S1 + S2

                with self.assertRaises(TypeError):
                    rdis3 = RaftDis.Create('rate', S1 + 2*S2)
                with self.assertRaises(TypeError):
                    rdis4 = RaftDis.Create(4, 'deps')
                with self.assertRaises(Exception):
                    rdis5 = RaftDis.Create(1, S1.o + 2* S2.i)

                stepsExo = rdis1._getStepsObjects()[0].__class__('rdis1b', rssys._getStepsObjects()[0], rdis1.AntiDependencies._GetStepsElems(), rdis1.K)

                rdis1b = RaftDis._FromStepsObject(stepsExo, self.mdl)
                self.assertEqual(rdis1b.K, rdis1.K)
                self.assertEqual(sorted([spec.getID() for spec in rdis1b.AntiDependencies._GetStepsElems()]), ['S1', 'S2', 'S2'])

            with self.assertRaises(Exception):
                rdis = RaftDis.Create(1, S1 + 2*S2)

    def testLinkSpeciesCreation(self):
        with self.mdl:
            lnkspec1 = LinkSpecies.Create(1)

            self.assertEqual(lnkspec1.Dcst, 1)

            with self.assertRaises(Exception):
                lnkspec1.Dcst = 2

            with self.assertRaises(TypeError):
                lnkspec3 = LinkSpecies.Create('dcst')

            stepsExo = lnkspec1._getStepsObjects()[0].__class__('lnkspec1b', self.mdl._getStepsObjects()[0], lnkspec1.Dcst)

            lnkspec1b = LinkSpecies._FromStepsObject(stepsExo, self.mdl)
            self.assertEqual(lnkspec1b.Dcst, lnkspec1.Dcst)

        with self.assertRaises(Exception):
            lnkspec = LinkSpecies.Create(1)

    def testVesSReacCreation(self):
        r = ReactionManager()

        with self.mdl:
            S1, S2, S3 = Species.Create()
            vssys = VesicleSurfaceSystem()

            with vssys:
                S1.v >r[1]> S2.v
                r[1].K = 1

            with vssys:
                S1.v >r[1]> S2.v
                r[1].K = 1
                r[1].Dependencies = S3.v

            with vssys:
                S1.v >r[1]> S2.v
                r[1].K = 1
                r[1].Dependencies = S3.v
                r[1].Immobilization = IMMOBILIZING

            # With patch species
            with vssys:
                S1.v + S3.s >r[1]> S2.v + S3.s
                r[1].K = 1
                r[1].Dependencies = S3.v
                r[1].Immobilization = IMMOBILIZING
                r[1].MaxDistance = 1e-6

            # Load from steps object
            so = r[1]._getStepsObjects()[0]
            so1 = so.__class__(
                'reac1b', vssys._getStepsObjects()[0], 
                so.getOLHS(), so.getSLHS(), so.getVLHS(), so.getLLHS(), 
                so.getLRHS(), so.getVRHS(), so.getSRHS(), so.getORHS(), so.getIRHS(), 
                so.getVDeps(), so.getKcst(), so.getImmobilization(), so.getMaxDistance()
            )
            with vssys:
                reac1b = Reaction._FromStepsObject(so1, self.mdl)
                self.assertEqual(set(reac1b.lhs), set(r[1].lhs))
                self.assertEqual(set(reac1b.rhs), set(r[1].rhs))
                self.assertEqual(reac1b.K, r[1].K)
                self.assertEqual(set(reac1b.Dependencies._toReactionSide()), set(r[1].Dependencies._toReactionSide()))
                self.assertEqual(reac1b.Immobilization, r[1].Immobilization)
                self.assertEqual(reac1b.MaxDistance, r[1].MaxDistance)

            # Everything on patch surface (strange but allowed)
            with vssys:
                S1.s + S2.s >r[1]> 2 * S1.s
                r[1].K = 1

            # Wrong parameter types
            for wdep in ['deps', 42]:
                with self.assertRaises(TypeError):
                    with vssys:
                        S1.v >r[1]> S1.v
                        r[1].K = 1
                        r[1].Dependencies = wdep

            for val, exc in [(True, TypeError), ('mob', TypeError), (-5, TypeError)]:
                with self.assertRaises(exc):
                    with vssys:
                        S1.v >r[1]> S1.v
                        r[1].K = 1
                        r[1].Immobilization = val

            for val, exc in [('42', TypeError), (-1e-6, ValueError)]:
                with self.assertRaises(exc):
                    with vssys:
                        S1.v + S2.s >r[1]> S1.v + S2.s
                        r[1].K = 1
                        r[1].MaxDistance = val

            # deps not on vesicle
            for wdep in [S3, S3.o, S3.r, S3.s, S3.i]:
                with self.assertRaises(Exception):
                    with vssys:
                        S1.v >r[1]> S1.v
                        r[1].K = 1
                        r[1].Dependencies = wdep

            # invalid location
            for val in [S2.i, S2.r, S2]:
                with self.assertRaises(Exception):
                    with vssys:
                        S1.v + val >r[1]> S1.v + S2.o
                        r[1].K = 1

            # Need to declare in VesSurfSys
            with self.assertRaises(Exception):
                S1.v >r[1]> S2.v
                S1.K = 1

            # TODO Reaction parameters as CompDepFunc

    def testRaftSReacCreation(self):
        r = ReactionManager()

        with self.mdl:
            S1, S2, S3 = Species.Create()
            rssys = RaftSurfaceSystem()

            with rssys:
                S1.r >r[1]> S2.r
                r[1].K = 1

            with rssys:
                S1.r >r[1]> S2.r
                r[1].K = 1
                r[1].Dependencies = S3.r

            with rssys:
                S1.r >r[1]> S2.r
                r[1].K = 1
                r[1].Dependencies = S3.r
                r[1].Immobilization = IMMOBILIZING

            # With patch species
            with rssys:
                S1.r + S3.s >r[1]> S2.r + S3.s
                r[1].K = 1
                r[1].Dependencies = S3.r
                r[1].Immobilization = IMMOBILIZING

            # Load from steps object
            so = r[1]._getStepsObjects()[0]
            so1 = so.__class__(
                'reac1b', rssys._getStepsObjects()[0], 
                so.getILHS(), so.getOLHS(), so.getSLHS(), so.getRLHS(), 
                so.getRRHS(), so.getSRHS(), so.getORHS(), so.getIRHS(), 
                so.getRDeps(), so.getAntiRDeps(), so.getKcst(), so.getImmobilization()
            )
            with rssys:
                reac1b = Reaction._FromStepsObject(so1, self.mdl)
                self.assertEqual(set(reac1b.lhs), set(r[1].lhs))
                self.assertEqual(set(reac1b.rhs), set(r[1].rhs))
                self.assertEqual(reac1b.K, r[1].K)
                self.assertEqual(set(reac1b.Dependencies._toReactionSide()), set(r[1].Dependencies._toReactionSide()))
                self.assertEqual(reac1b.Immobilization, r[1].Immobilization)

            # Everything on patch surface (strange but allowed)
            with rssys:
                S1.s + S2.s >r[1]> 2 * S1.s
                r[1].K = 1

            # Wrong parameter types
            for wdep in ['deps', 42]:
                with self.assertRaises(TypeError):
                    with rssys:
                        S1.r >r[1]> S1.r
                        r[1].K = 1
                        r[1].Dependencies = wdep

            for val, exc in [(True, TypeError), ('mob', TypeError), (-5, TypeError)]:
                with self.assertRaises(exc):
                    with rssys:
                        S1.r >r[1]> S1.r
                        r[1].K = 1
                        r[1].Immobilization = val

            # deps not on raft
            for wdep in [S3, S3.o, S3.v, S3.s, S3.i]:
                with self.assertRaises(Exception):
                    with rssys:
                        S1.r >r[1]> S1.r
                        r[1].K = 1
                        r[1].Dependencies = wdep

            # invalid location
            for val in [S2.v, S2]:
                with self.assertRaises(Exception):
                    with rssys:
                        S1.r + val >r[1]> S1.r + S2.o
                        r[1].K = 1

            # Need to declare in RaftSurfSys
            with self.assertRaises(Exception):
                S1.r >r[1]> S2.r
                S1.K = 1

            # TODO Reaction parameters as CompDepFunc

    def testVesSDiffCreation(self):

        with self.mdl:
            S1, S2, S3 = Species.Create()
            vssys = VesicleSurfaceSystem()

            with vssys:
                diff = Diffusion(S2, 1)
                self.assertEqual(diff.Dcst, 1)

                diff1 = Diffusion.Create(S1, 1)
                self.assertEqual(diff1.Dcst, 1)

            diff1.Dcst = 2
            self.assertEqual(diff1.Dcst, 2)

            # Load from steps object
            so = diff1._getStepsObjects()[0]
            so1 = so.__class__(
                'diff1b', vssys._getStepsObjects()[0], so.getLig(), so.getDcst()
            )
            with vssys:
                diff1b = Diffusion._FromStepsObject(so1, self.mdl)
                self.assertEqual(diff1b._elem, diff1._elem)
                self.assertEqual(diff1b.Dcst, diff1.Dcst)

            # Wrong parameter types
            for val in ['lig', 42]:
                with self.assertRaises(TypeError):
                    with vssys:
                        Diffusion(val, 1)

            for val in [S1, '42']:
                with self.assertRaises(TypeError):
                    with vssys:
                        Diffusion(S1, val)

            # Need to declare in VesSurfSys
            with self.assertRaises(Exception):
                Diffusion(S1, 1)

            # TODO Diffusion parameters as CompDepFunc

    def testVesicleBindCreation(self):
        r = ReactionManager()

        with self.mdl:
            S1, S2, S3 = Species.Create()
            L1, L2, L3 = LinkSpecies.Create()
            vsys = VolumeSystem()

            ves1, ves2 = Vesicle.Create(1, 2) 
            raft1 = Raft.Create(1)

            with vsys:
                vb1 = VesicleBind.Create((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7)
                self.assertEqual(vb1.LengthMin, 4e-7)
                self.assertEqual(vb1.LengthMax, 5e-7)
                self.assertEqual(vb1.K, 0)
                self.assertEqual(vb1.Immobilization, NO_EFFECT)
                self.assertEqual(vb1.Dependencies[0], None)
                self.assertEqual(vb1.Dependencies[1], None)

            with vsys:
                vb2 = VesicleBind.Create((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                self.assertEqual(vb2.LengthMin, 4e-7)
                self.assertEqual(vb2.LengthMax, 5e-7)
                self.assertEqual(vb2.K, 1e-12)
                self.assertEqual(vb2.Immobilization, IMMOBILIZING)
                self.assertEqual(set(vb2.Dependencies[0]), set(S3.v + L3.v))
                self.assertEqual(set(vb2.Dependencies[1]), set(S2.v + L2.v))

                vb2.LengthMin = 6e-7
                vb2.LengthMax = 7e-7
                vb2.K = 2e-12
                vb2.Immobilization = NO_EFFECT
                vb2.Dependencies = (S1.v + L1.v, S3.v + L3.v)

            self.assertEqual(vb2.LengthMin, 6e-7)
            self.assertEqual(vb2.LengthMax, 7e-7)
            self.assertEqual(vb2.K, 2e-12)
            self.assertEqual(vb2.Immobilization, NO_EFFECT)
            self.assertEqual(set(vb2.Dependencies[0]), set(S1.v + L1.v))
            self.assertEqual(set(vb2.Dependencies[1]), set(S3.v + L3.v))

            # Cant change parameters after creation
            with self.assertRaises(Exception):
                vb1.K = 5
            with self.assertRaises(Exception):
                vb1.LengthMin = 5
            with self.assertRaises(Exception):
                vb1.LengthMax = 5
            with self.assertRaises(Exception):
                vb1.Immobilization = IMMOBILIZING
            with self.assertRaises(Exception):
                vb1.Dependencies = (None, None)

            # Wrong parameter types
            with vsys:
                with self.assertRaises(TypeError):
                    VesicleBind('ves1 ves2', (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind(('ves1', 'ves2'), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((raft1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), ('S1', 'S2'), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (L1, L2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), '(S1, S2)', (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (S1, S2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), '(L1, L2)', 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), ('L1', 'L2'), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (L1, L2), '4e-7', 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (L1, L2), 4e-7, '5e-7', 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, '1e-12', IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, 1, (S3.v + L3.v, S2.v + L2.v))
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, '(S3.v + L3.v, S2.v + L2.v)')
                with self.assertRaises(TypeError):
                    VesicleBind((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, ('S3.v + L3.v', 'S2.v + L2.v'))

                vb3 = VesicleBind.Create((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.v + L3.v, S2.v + L2.v))

                with self.assertRaises(TypeError):
                    vb3.LengthMin = '6e-7'
                with self.assertRaises(TypeError):
                    vb3.LengthMax = '7e-7'
                with self.assertRaises(TypeError):
                    vb3.K = '2e-12'
                with self.assertRaises(TypeError):
                    vb3.Immobilization = '-1'
                with self.assertRaises(TypeError):
                    vb3.Immobilization = -2
                with self.assertRaises(TypeError):
                    vb3.Dependencies = '(S1.v + L1.v, S3.v + L3.v)'
                with self.assertRaises(TypeError):
                    vb3.Dependencies = S1.v + L1.v

            # Wrong position for dependencies
            with self.assertRaises(Exception):
                with vsys:
                    vb4 = VesicleBind.Create((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.i + L3.v, S2.o + L2.v))


            # Load from steps object
            so = vb3._getStepsObjects()[0]
            so1 = so.__class__(
                'vb3b', vsys._getStepsObjects()[0], 
                so.getVesicle1(), so.getReactant1(), so.getVesicle2(), so.getReactant2(), 
                so.getProduct1(), so.getProduct2(), so.getLengthMax(), so.getLengthMin(), so.getVDeps1(), 
                so.getVDeps2(), so.getLDeps1(), so.getLDeps2(), so.getKcst(), so.getImmobilization()
            )
            with vsys:
                vb3b = VesicleBind._FromStepsObject(so1, self.mdl)
                self.assertEqual(vb3b.Vesicles, vb3.Vesicles)
                self.assertEqual(vb3b._reactants, vb3._reactants)
                self.assertEqual(vb3b._linkProducts, vb3._linkProducts)
                self.assertEqual(vb3b.LengthMin, vb3.LengthMin)
                self.assertEqual(vb3b.LengthMax, vb3.LengthMax)
                self.assertEqual(vb3b.K, vb3.K)
                self.assertEqual(vb3b.Immobilization, vb3.Immobilization)
                self.assertEqual(set(vb3b.Dependencies[0]), set(vb3.Dependencies[0]))
                self.assertEqual(set(vb3b.Dependencies[1]), set(vb3.Dependencies[1]))

            # Need to declare in VesSurfSys
            with self.assertRaises(Exception):
                vb5 = VesicleBind.Create((ves1, ves2), (S1, S2), (L1, L2), 4e-7, 5e-7, 1e-12, IMMOBILIZING, (S3.i + L3.v, S2.o + L2.v))

    def testVesicleUnbindCreation(self):
        r = ReactionManager()

        with self.mdl:
            S1, S2, S3 = Species.Create()
            L1, L2, L3 = LinkSpecies.Create()
            vsys = VolumeSystem()

            ves1, ves2 = Vesicle.Create(1, 2) 
            raft1 = Raft.Create(1)

            with vsys:
                vb1 = VesicleUnbind.Create((ves1, ves2), (L1, L2), (S1, S2))
                self.assertEqual(vb1.K, 0)
                self.assertEqual(vb1.Immobilization, NO_EFFECT)

            with vsys:
                vb2 = VesicleUnbind.Create((ves1, ves2), (L1, L2), (S1, S2), 5, IMMOBILIZING)
                self.assertEqual(vb2.K, 5)
                self.assertEqual(vb2.Immobilization, IMMOBILIZING)

                vb2.K = 2e-12
                vb2.Immobilization = MOBILIZING

            self.assertEqual(vb2.K, 2e-12)
            self.assertEqual(vb2.Immobilization, MOBILIZING)

            # Cant change parameters after creation
            with self.assertRaises(Exception):
                vb1.K = 5
            with self.assertRaises(Exception):
                vb1.Immobilization = IMMOBILIZING

            # Wrong parameter types
            with vsys:
                with self.assertRaises(TypeError):
                    VesicleUnbind('(ves1, ves2)', (L1, L2), (S1, S2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind(('ves1', 'ves2'), (L1, L2), (S1, S2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((raft1, ves2), (L1, L2), (S1, S2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), '(L1, L2)', (S1, S2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), ('L1', 'L2'), (S1, S2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (S1, L2), (S1, S2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (L1, L2), '(S1, S2)', 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (L1, L2), ('S1', 'S2'), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (L1, L2), (S1, L2), 5, MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (L1, L2), (S1, S2), '5', MOBILIZING)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (L1, L2), (S1, S2), 5, 1)
                with self.assertRaises(TypeError):
                    VesicleUnbind((ves1, ves2), (L1, L2), (S1, S2), 5, -5)

                vb3 = VesicleUnbind.Create((ves1, ves2), (L1, L2), (S1, S2), 5, MOBILIZING)

                with self.assertRaises(TypeError):
                    vb3.K = '2e-12'
                with self.assertRaises(TypeError):
                    vb3.Immobilization = '-1'
                with self.assertRaises(TypeError):
                    vb3.Immobilization = -2

            # Load from steps object
            so = vb3._getStepsObjects()[0]
            so1 = so.__class__(
                'vb3b', vsys._getStepsObjects()[0], 
                so.getLink1(), so.getLink2(),
                so.getVesicle1(), so.getProduct1(), so.getVesicle2(), so.getProduct2(), 
                so.getKcst(), so.getImmobilization()
            )
            with vsys:
                vb3b = VesicleUnbind._FromStepsObject(so1, self.mdl)
                self.assertEqual(vb3b.Vesicles, vb3.Vesicles)
                self.assertEqual(vb3b._linkReactants, vb3._linkReactants)
                self.assertEqual(vb3b._products, vb3._products)
                self.assertEqual(vb3b.K, vb3.K)
                self.assertEqual(vb3b.Immobilization, vb3.Immobilization)

            # VesicleUnbind cannot be IMMOBILIZING
            with self.assertRaises(Exception):
                with vsys:
                    vb3c = VesicleUnbind.Create((ves1, ves2), (L1, L2), (S1, S2), 5, IMMOBILIZING)

            # Need to declare in VesSurfSys
            with self.assertRaises(Exception):
                vb4 = VesicleUnbind.Create((ves1, ves2), (L1, L2), (S1, S2), 5, MOBILIZING)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(VesicleObjectsDeclaration, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


