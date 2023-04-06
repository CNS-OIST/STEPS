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

""" Default model to be used for testing """

import os
import sys
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *
from steps.sim import SimPathInvalidPath, SimPathSolverMissingMethod, SolverCallError

import steps.API_1.model as smodel
import steps.API_1.geom as sgeom
import steps.API_1.rng as srng
import steps.API_1.solver as ssolver
import steps.API_1.utilities.meshio as smeshio
import steps.API_1.utilities.meshctrl as smeshctrl

FILEDIR = os.path.dirname(os.path.abspath(__file__))

Avogad = 6.02214076e23


class TestModelFramework(unittest.TestCase):

    def setUp(self):
        self.seed = 12345

        self.endTime = 1
        self.shortEndTime = 0.1
        self.deltaT = 0.01
        self.nbRuns = 3

        self.tolerance = 0.15
        self.countRefVal = 150

        # Mdl parameters

        self.r01f = 2000
        self.r02f = 5000
        self.r03f = 4000000
        self.r04f = 5100000
        self.r05f = 4200000
        self.r06f = 5300000
        self.r07f = 10
        self.cr01f = 2451

        self.r01b = 100
        self.r02b = 2500
        self.r03b = 1.0
        self.r04b = 1.5
        self.r05b = 1.1
        self.r06b = 1.4
        self.cr01b = 152

        self.initC1S1 = 500
        self.initC1S2 = 900
        self.initC1S3 = 1000
        self.initC2S3 = 10
        self.initC2S1 = 1000
        self.initC2S2 = 500
        self.initP1Ex = 1000
        self.initP1ExS1 = 11
        self.initP1ExS2 = 15
        self.initP1ExS1S2 = 17

        self.initC1CCsus11 = 30
        self.initC1CCsus12 = 15
        self.initC1CCsus22 = 60

        # Geom parameters 

        self.v1 = 1e-20
        self.v2 = 2e-20
        self.a1 = 3e-14

        self.useDist = False

    def assertSameData(self, oldDat, newDat, refVal):
        self.assertEqual(len(oldDat), len(newDat))
        for orun, nrun in zip(oldDat, newDat):
            self.assertEqual(len(orun), len(nrun))
            for orow, nrow in zip(orun, nrun):
                self.assertEqual(len(orow), len(nrow))
                for o, n, r in zip(orow, nrow, refVal):
                    if abs(o + n)/2 > r:
                        # Uncomment for displaying the difference between datasets
                        # import pylab
                        # if MPI._shouldWrite and abs(o - n) / ((o+n)/2)>= self.tolerance:
                            # print(o, n, abs(o - n) / ((o+n)/2))
                            # pylab.plot(range(orun.shape[0]), orun, 'b', range(nrun.shape[0]), nrun, 'r')
                            # pylab.show()
                        self.assertLess(abs(o - n) / ((o+n)/2), self.tolerance)


    def get_API2_Mdl(self, sas=True):
        # New interface code
        nmdl = Model.Create()
        r = ReactionManager()
        with nmdl:
            S1, S2, S3, Ex, ExS1, ExS2, ExS1S2 = Species.Create()

            sus1, sus2 = SubUnitState.Create()
            su = SubUnit.Create([sus1, sus2])
            CC = Complex.Create([su, su], statesAsSpecies=sas)
            self.sus1, self.sus2 = sus1, sus2

            vsys1, vsys2 = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

            with vsys1:
                S1 + S2 <r['vs1R1']> 2*S1
                r['vs1R1'].K = self.r01f, self.r01b

                with CC[...]:
                    sus1 + S1 <r['vs1CR1']> sus2
                    r['vs1CR1'].K = self.cr01f, self.cr01b

            with vsys2:
                S1 + S2 <r['vs2R2']> 2*S2
                r['vs2R2'].K = self.r02f, self.r02b

                # TODO Temp to prevent the DistTetOpSplit bug not supporting volume systems with different
                # numbers of species
                with CC[...]:
                    sus1 + S1 <r['vs1CR2']> sus2
                    r['vs1CR2'].K = self.cr01f, self.cr01b

            with ssys:
                (S1.i + Ex.s <r['ss1R3']> ExS1.s) + S2.o <r['ss1R4']> ExS1S2.s
                r['ss1R3'].K = self.r03f, self.r03b
                r['ss1R4'].K = self.r04f, self.r04b
                (S2.o + Ex.s <r['ss1R5']> ExS2.s) + S1.i <r['ss1R6']> ExS1S2.s
                r['ss1R5'].K = self.r05f, self.r05b
                r['ss1R6'].K = self.r06f, self.r06b

                ExS1S2.s >r['ss1R7']> S1.o + S2.i + Ex.s
                r['ss1R7'].K = self.r07f

        return nmdl

    def get_API1_Mdl(self):
        # Old interface code
        omdl = smodel.Model()

        S1 = smodel.Spec('S1', omdl)
        S2 = smodel.Spec('S2', omdl)
        S3 = smodel.Spec('S3', omdl)
        CCsus11 = smodel.Spec('CCsus11', omdl)
        CCsus12 = smodel.Spec('CCsus12', omdl)
        CCsus22 = smodel.Spec('CCsus22', omdl)
        Ex = smodel.Spec('Ex', omdl)
        ExS1 = smodel.Spec('ExS1', omdl)
        ExS2 = smodel.Spec('ExS2', omdl)
        ExS1S2 = smodel.Spec('ExS1S2', omdl)

        vsys1 = smodel.Volsys('vsys1', omdl)
        vsys2 = smodel.Volsys('vsys2', omdl)
        ssys = smodel.Surfsys('ssys', omdl)

        smodel.Reac('R01f', vsys1, lhs=[S1, S2], rhs=[S1, S1], kcst=self.r01f)
        smodel.Reac('R01b', vsys1, lhs=[S1, S1], rhs=[S1, S2], kcst=self.r01b)

        smodel.Reac('R02f', vsys2, lhs=[S1, S2], rhs=[S2, S2], kcst=self.r02f)
        smodel.Reac('R02b', vsys2, lhs=[S2, S2], rhs=[S1, S2], kcst=self.r02b)

        smodel.SReac('R03f', ssys, ilhs=[S1], slhs=[Ex], srhs=[ExS1], kcst=self.r03f)
        smodel.SReac('R03b', ssys, irhs=[S1], srhs=[Ex], slhs=[ExS1], kcst=self.r03b)

        smodel.SReac('R04f', ssys, olhs=[S2], slhs=[ExS1], srhs=[ExS1S2], kcst=self.r04f)
        smodel.SReac('R04b', ssys, orhs=[S2], srhs=[ExS1], slhs=[ExS1S2], kcst=self.r04b)

        smodel.SReac('R05f', ssys, olhs=[S2], slhs=[Ex], srhs=[ExS2], kcst=self.r05f)
        smodel.SReac('R05b', ssys, orhs=[S2], srhs=[Ex], slhs=[ExS2], kcst=self.r05b)

        smodel.SReac('R06f', ssys, ilhs=[S1], slhs=[ExS2], srhs=[ExS1S2], kcst=self.r06f)
        smodel.SReac('R06b', ssys, irhs=[S1], srhs=[ExS2], slhs=[ExS1S2], kcst=self.r06b)

        smodel.SReac('R07f', ssys, slhs=[ExS1S2], srhs=[Ex], irhs=[S2], orhs=[S1], kcst=self.r07f)

        smodel.Reac('CR01_1f', vsys1, lhs=[S1, CCsus11], rhs=[CCsus12], kcst=2*self.cr01f)
        smodel.Reac('CR01_1b', vsys1, lhs=[CCsus12], rhs=[S1, CCsus11], kcst=self.cr01b)
        smodel.Reac('CR01_2f', vsys1, lhs=[S1, CCsus12], rhs=[CCsus22], kcst=1*self.cr01f)
        smodel.Reac('CR01_2b', vsys1, lhs=[CCsus22], rhs=[S1, CCsus12], kcst=2*self.cr01b)

        return omdl

    def get_API2_Geom(self, nmdl):
        # New interface code
        ngeom = Geometry()
        with ngeom:
            comp1, comp2 = Compartment.Create(
                Params(nmdl.vsys1, self.v1), 
                Params(nmdl.vsys2, self.v2)
            )
            patch = Patch.Create(comp1, comp2, nmdl.ssys, self.a1)

        return ngeom

    def get_API1_Geom(self, omdl):
        # Old interface code
        ogeom = sgeom.Geom()

        comp1 = sgeom.Comp('comp1', ogeom)
        comp1.setVol(self.v1)
        comp1.addVolsys('vsys1')

        comp2 = sgeom.Comp('comp2', ogeom)
        comp2.setVol(self.v2)
        comp2.addVolsys('vsys2')

        patch = sgeom.Patch('patch', ogeom, comp1, comp2)
        patch.addSurfsys('ssys')
        patch.setArea(self.a1)

        return ogeom

    def init_API2_sim(self, sim):
        sim.comp1.S1.Count = self.initC1S1
        sim.comp1.S2.Count = self.initC1S2
        sim.comp1.CC[self.sus1, self.sus1].Count = self.initC1CCsus11
        sim.comp1.CC[self.sus1, self.sus2].Count = self.initC1CCsus12
        sim.comp1.CC[self.sus2, self.sus2].Count = self.initC1CCsus22
        sim.comp2.S1.Count = self.initC2S1
        sim.comp2.S2.Count = self.initC2S2
        sim.patch.Ex.Count = self.initP1Ex
        sim.patch.ExS1.Count = self.initP1ExS1
        sim.patch.ExS2.Count = self.initP1ExS2
        sim.patch.ExS1S2.Count = self.initP1ExS1S2

    def init_API1_sim(self, sim):
        sim.setCompCount('comp1', 'S1', self.initC1S1)
        sim.setCompCount('comp1', 'S2', self.initC1S2)
        sim.setCompCount('comp1', 'CCsus11', self.initC1CCsus11)
        sim.setCompCount('comp1', 'CCsus12', self.initC1CCsus12)
        sim.setCompCount('comp1', 'CCsus22', self.initC1CCsus22)
        sim.setCompCount('comp2', 'S1', self.initC2S1)
        sim.setCompCount('comp2', 'S2', self.initC2S2)
        sim.setPatchCount('patch', 'Ex', self.initP1Ex)
        sim.setPatchCount('patch', 'ExS1', self.initP1ExS1)
        sim.setPatchCount('patch', 'ExS2', self.initP1ExS2)
        sim.setPatchCount('patch', 'ExS1S2', self.initP1ExS1S2)

    def assertAlmostEqualWThresh(self, v1, v2, thresh):
        mean = (v1 + v2) / 2
        self.assertAlmostEqual(v1, v2, delta=abs(thresh * mean))

    def _test_API2_SimPathGetSyntax(self, sim, init=True, usingMesh=False, extents=None):
        if init:
            self.init_API2_sim(sim)

        sus1, sus2 = self.sus1, self.sus2

        if not self.useDist:
            self.assertEqual(sim.comp1.Vol, self.v1)
            self.assertEqual(sim.comp2.Vol, self.v2)
            self.assertEqual(sim.patch.Area, self.a1)

        # Solver global values
        if not self.useDist:
            self.assertTrue(sim.A0 > 0)
        if usingMesh:
            self.assertEqual(sim.EfieldDT, self.efielddt)
            self.assertEqual(sim.Temp, self.solvTemp)

        # Species
        specData1 = [
            (sim.comp1.S1, self.initC1S1, self.v1),
            (sim.comp1.S2, self.initC1S2, self.v1),
            (sim.comp1.CC[sus1, sus1], self.initC1CCsus11, self.v1),
            (sim.comp1.CC[sus1, sus2], self.initC1CCsus12, self.v1),
            (sim.comp1.CC[sus2, sus2], self.initC1CCsus22, self.v1),
            (sim.comp1.CC.sus2, self.initC1CCsus12 + 2*self.initC1CCsus22, self.v1),
            (sim.comp2.S1, self.initC2S1, self.v2),
            (sim.comp2.S2, self.initC2S2, self.v2),
        ]
        for path, cnt, vol in specData1:
            with self.subTest(path=path, cnt=cnt, vol=vol):
                # TODO Remove condition once vesicle changes are merged
                pathCnt = path.Count
                pathCnc = path.Conc
                if MPI.rank == 0 or not self.useDist:
                    self.assertEqual(pathCnt, cnt)
                    self.assertAlmostEqualWThresh(pathCnc, cnt / Avogad / (1e3 * vol), self.tolerance)
                if not self.useDist:
                    self.assertAlmostEqualWThresh(path.Amount, cnt / Avogad, self.tolerance)
                    self.assertFalse(path.Clamped)

        specData2 = [
            (sim.patch.Ex, self.initP1Ex),
            (sim.patch.ExS1, self.initP1ExS1),
            (sim.patch.ExS2, self.initP1ExS2),
            (sim.patch.ExS1S2, self.initP1ExS1S2),
        ]
        for path, cnt in specData2:
            with self.subTest(path=path, cnt=cnt):
                self.assertEqual(path.Count, cnt)
                if not self.useDist:
                    self.assertAlmostEqualWThresh(path.Amount, cnt / Avogad, self.tolerance)
                    self.assertFalse(path.Clamped)

        # Reactions
        reacData = [
            (sim.comp1.vs1R1, 'fwd', self.r01f, self.initC1S1 * self.initC1S2, 1e3 * self.v1 * Avogad),
            (sim.comp1.vs1R1, 'bkw', self.r01b, self.initC1S1 * (self.initC1S1 - 1), 1e3 * self.v1 * Avogad),
            (sim.comp2.vs2R2, 'fwd', self.r02f, self.initC2S1 * self.initC2S2, 1e3 * self.v2 * Avogad),
            (sim.comp2.vs2R2, 'bkw', self.r02b, self.initC2S2 * (self.initC2S2 - 1), 1e3 * self.v2 * Avogad),
            (sim.patch.ss1R3, 'fwd', self.r03f, self.initC1S1 * self.initP1Ex, 1e3 * self.v1 * Avogad),
            (sim.patch.ss1R3, 'bkw', self.r03b, self.initP1ExS1, 1),
            (sim.patch.ss1R4, 'fwd', self.r04f, self.initP1ExS1 * self.initC2S2, 1e3 * self.v2 * Avogad),
            (sim.patch.ss1R4, 'bkw', self.r04b, self.initP1ExS1S2, 1),
            (sim.patch.ss1R5, 'fwd', self.r05f, self.initC2S2 * self.initP1Ex, 1e3 * self.v2 * Avogad),
            (sim.patch.ss1R5, 'bkw', self.r05b, self.initP1ExS2, 1),
            (sim.patch.ss1R6, 'fwd', self.r06f, self.initP1ExS2 * self.initC1S1, 1e3 * self.v1 * Avogad),
            (sim.patch.ss1R6, 'bkw', self.r06b, self.initP1ExS1S2, 1),
            (sim.patch.ss1R7, None, self.r07f, self.initP1ExS1S2, 1),
        ]
        if not self.useDist:
            if extents is None:
                extents = [0]*len(reacData)
            for (path, specif, k, h, vs), extent in zip(reacData, extents):
                with self.subTest(path=path, specif=specif, k=k, h=h, vs=vs):
                    if specif is not None:
                        path = path[specif]
                    self.assertEqual(path.K, k)
                    self.assertTrue(path.Active)
                    self.assertEqual(path.Extent, extent)
                    if not usingMesh:
                        self.assertAlmostEqualWThresh(path.H, h, self.tolerance)
                        self.assertAlmostEqualWThresh(path.C, k / vs, self.tolerance)
                        self.assertAlmostEqualWThresh(path.A, h * k / vs, self.tolerance)
            self.assertEqual(set(sim.comp1.vs1CR1['fwd'].K), set([self.cr01f, 2 * self.cr01f]))
            self.assertEqual(set(sim.comp1.vs1CR1['bkw'].K), set([self.cr01b, 2 * self.cr01b]))
        if usingMesh:
            if not self.useDist:
                self.assertTrue(sim.patch.ss1RVDep01['fwd'].Active)
            self.assertFalse(sim.diffb.S1.DiffusionActive)
            if not self.useDist:
                self.assertFalse(sim.sdiffb.Ex.DiffusionActive)

        # Special methods
        self.assertEqual(
            sum(sim.patch.MATCH('Ex*').Count), 
            self.initP1Ex + self.initP1ExS1 + self.initP1ExS2 + self.initP1ExS1S2
        )
        if usingMesh:
            tet = self.newGeom.comp1.tets[0]
            if not self.useDist:
                self.assertEqual(sim.TET(tet).Vol, tet.Vol)
                self.assertAlmostEqualWThresh(sim.TET(tet).S1.Conc, sim.TET(tet).S1.Amount / tet.Vol * 1e-3, self.tolerance)
                self.assertEqual(sim.TET(tet).S1.Clamped, False)

                tri = self.newGeom.patch.tris[0]
                self.assertEqual(sim.TRI(tri).Area, tri.Area)

                self.assertEqual(sim.comp1.Diffc1S1.Active, True)

                self.assertEqual(sim.comp1.Diffc1S1.D, self.c1DS1)
                self.assertEqual(sim.comp1.Diffc1S2.D, self.c1DS2)
                self.assertEqual(sim.comp2.Diffc2S1.D, self.c2DS1)
                self.assertEqual(sim.comp2.Diffc2S2.D, self.c2DS2)

                CC, sus1, sus2, = self.newMdl.CC, self.newMdl.sus1, self.newMdl.sus2
                self.assertEqual(set(sim.comp1.Diffc1CC.D), set([self.c1DCC(state) for state in self.newMdl.CC]))
                self.assertEqual(set(sim.comp1.Diffc1CC[sus2, :].D), set([self.c1DCC(state) for state in CC[sus2,:]]))
                self.assertEqual(sim.comp1.Diffc1CC[self.newMdl.sus1, self.newMdl.sus1].D, self.c1DCC(CC[sus1, sus1]))
                self.assertEqual(sim.comp1.Diffc1CC[self.newMdl.sus1, self.newMdl.sus2].D, self.c1DCC(CC[sus1, sus2]))
                self.assertEqual(sim.comp1.Diffc1CC[self.newMdl.sus2, self.newMdl.sus2].D, self.c1DCC(CC[sus2, sus2]))

            self.assertAlmostEqualWThresh(
                sum(sim.TETS(self.newGeom.comp1.tets).S1.Count),
                sim.comp1.S1.Count,
                self.tolerance
            )
            if not self.useDist:
                self.assertAlmostEqualWThresh(
                    sum(sim.TETS(self.newGeom.comp1.tets).S1.Amount),
                    sim.comp1.S1.Amount,
                    self.tolerance
                )
            self.assertAlmostEqualWThresh(
                sum(sim.TRIS(self.newGeom.patch.tris).ExS1.Count),
                sim.patch.ExS1.Count,
                self.tolerance
            )
            self.assertAlmostEqualWThresh(
                sum(sim.TET(t).S1.Count for t in self.newGeom.comp1.tets),
                sim.comp1.S1.Count,
                self.tolerance
            )
            self.assertAlmostEqualWThresh(
                sum(sim.TRI(t).ExS1.Count for t in self.newGeom.patch.tris),
                sim.patch.ExS1.Count,
                self.tolerance
            )
            # TODO Remove the condition after changes from vesicle branch have been merged
            tetVal = sum(sim.TETS().S1.Count)
            compVal = sim.comp1.S1.Count + sim.comp2.S1.Count
            if MPI._shouldWrite or not self.useDist:
                self.assertAlmostEqualWThresh(
                    tetVal,
                    compVal,
                    self.tolerance
                )
            if not self.useDist:
                self.assertAlmostEqualWThresh(
                    sum(sim.TRIS(self.newGeom.patch.tris).Area),
                    sim.patch.Area,
                    self.tolerance
                )
            for v in sim.VERTS(self.newGeom.membrane.tris[0].verts).V:
                self.assertAlmostEqual(v, self.membPot)
            self.assertAlmostEqual(
                sim.VERT(self.newGeom.membrane.tris[0].verts[0]).V,
                self.membPot
            )
            if not self.useDist:
                self.assertFalse(sim.VERT(self.newGeom.membrane.tris[0].verts[0]).VClamped)

            sim.TRIS()
            sim.VERTS()

            with self.assertRaises(SimPathInvalidPath):
                sim.comp1.TETS().S1.Count
            with self.assertRaises(SimPathInvalidPath):
                sim.comp1.TET(0).S1.Count
            with self.assertRaises(SimPathInvalidPath):
                sim.patch.TRIS().ExS1.Count
            with self.assertRaises(SimPathInvalidPath):
                sim.patch.TRI(0).ExS1.Count
            with self.assertRaises(SimPathInvalidPath):
                sim.membrane.VERTS().V
            with self.assertRaises(SimPathInvalidPath):
                sim.membrane.VERT(0).V


            tri1 = self.newGeom.patch.tris[0]
            tet1 = [tet for tet in tri1.tetNeighbs if tet in self.newGeom.comp1.tets][0]
            tet2 = [tet for tet in tri1.tetNeighbs if tet in self.newGeom.comp2.tets][0]

            self.assertEqual(sim.TET(tet1).V, self.membPot)
            if not self.useDist:
                self.assertEqual(sim.TET(tet1).VClamped, False)
            self.assertAlmostEqual(sim.TRI(tri1).V, self.membPot)
            if not self.useDist:
                self.assertEqual(sim.TRI(tri1).VClamped, False)

            if not self.useDist:
                # Species
                cnttmp = sim.TET(tet1).S1.Count
                self.assertAlmostEqualWThresh(sim.TET(tet1).S1.Amount, cnttmp / Avogad, self.tolerance)
                self.assertAlmostEqualWThresh(sim.TET(tet1).S1.Conc, cnttmp / Avogad / (1e3 * tet1.Vol), self.tolerance)
                self.assertFalse(sim.TET(tet1).S1.Clamped)

                cnttmp = sim.TRI(tri1).Ex.Count
                self.assertAlmostEqualWThresh(sim.TRI(tri1).Ex.Amount, cnttmp / Avogad, self.tolerance)
                self.assertFalse(sim.TRI(tri1).Ex.Clamped)

                # Diffusion
                self.assertEqual(sim.TET(tet1).Diffc1S1.D, self.c1DS1)
                self.assertEqual(sim.TET(tet1).Diffc1S1.Active, True)
                sim.TET(tet1).Diffc1S1.A

                self.assertEqual(sim.TRI(tri1).DiffpatchEx.D, self.patchDEx)

                # Reactions
                reacData = [
                    (sim.TET(tet1).vs1R1, 'fwd', self.r01f, sim.TET(tet1).S1.Count * sim.TET(tet1).S2.Count, 1e3 * tet1.Vol * Avogad),
                    (sim.TET(tet1).vs1R1, 'bkw', self.r01b, sim.TET(tet1).S1.Count * (sim.TET(tet1).S1.Count - 1), 1e3 * tet1.Vol * Avogad),
                    (sim.TET(tet2).vs2R2, 'fwd', self.r02f, sim.TET(tet2).S1.Count * sim.TET(tet2).S2.Count, 1e3 * tet2.Vol * Avogad),
                    (sim.TET(tet2).vs2R2, 'bkw', self.r02b, sim.TET(tet2).S2.Count * (sim.TET(tet2).S2.Count - 1), 1e3 * tet2.Vol * Avogad),
                    (sim.TRI(tri1).ss1R3, 'fwd', self.r03f, sim.TET(tet1).S1.Count * sim.TRI(tri1).Ex.Count, 1e3 * tet1.Vol * Avogad),
                    (sim.TRI(tri1).ss1R3, 'bkw', self.r03b, sim.TRI(tri1).ExS1.Count, 1),
                    (sim.TRI(tri1).ss1R4, 'fwd', self.r04f, sim.TRI(tri1).ExS1.Count * sim.TET(tet2).S2.Count, 1e3 * tet2.Vol * Avogad),
                    (sim.TRI(tri1).ss1R4, 'bkw', self.r04b, sim.TRI(tri1).ExS1S2.Count, 1),
                    (sim.TRI(tri1).ss1R5, 'fwd', self.r05f, sim.TET(tet2).S2.Count * sim.TRI(tri1).Ex.Count, 1e3 * tet2.Vol * Avogad),
                    (sim.TRI(tri1).ss1R5, 'bkw', self.r05b, sim.TRI(tri1).ExS2.Count, 1),
                    (sim.TRI(tri1).ss1R6, 'fwd', self.r06f, sim.TRI(tri1).ExS2.Count * sim.TET(tet1).S1.Count, 1e3 * tet1.Vol * Avogad),
                    (sim.TRI(tri1).ss1R6, 'bkw', self.r06b, sim.TRI(tri1).ExS1S2.Count, 1),
                    (sim.TRI(tri1).ss1R7, None, self.r07f, sim.TRI(tri1).ExS1S2.Count, 1),
                ]
                for path, specif, k, h, vs in reacData:
                    with self.subTest(path=path, specif=specif, k=k, h=h, vs=vs):
                        if specif is not None:
                            path = path[specif]
                        self.assertEqual(path.K, k)
                        self.assertTrue(path.Active)
                        # Extent not implemented in STEPS for tets and tris
                        # self.assertEqual(path.Extent, 0)
                        self.assertAlmostEqualWThresh(path.H, h, self.tolerance)
                        self.assertAlmostEqualWThresh(path.C, k / vs, self.tolerance)
                        self.assertAlmostEqualWThresh(path.A, h * k / vs, self.tolerance)

                # Currents
                ohmcurr = sim.TRI(tri1).Chan1_Ohm_I.I
                ghkcurr = sim.TRI(tri1).Chan1_GHK_I.I
                self.assertAlmostEqualWThresh(sim.TRI(tri1).I, ohmcurr + ghkcurr, self.tolerance)

                # VDepSReac
                self.assertTrue(sim.patch.ss1RVDep01['fwd'].Active)
                self.assertTrue(sim.TRI(tri1).ss1RVDep01['fwd'].Active)

                #ROI
                self.assertEqual(sim.comp1ROI.Vol, sim.comp1.Vol)
                self.assertEqual(sim.patchROI.Area, sim.patch.Area)

                self.assertEqual(sim.comp1ROI.S1.Count, sim.comp1.S1.Count)
                self.assertEqual(sim.comp1ROI.S1.Amount, sim.comp1.S1.Amount)
                self.assertEqual(sim.comp1ROI.S1.Conc, sim.comp1.S1.Conc)
                if sim._solverStr != 'TetOpSplit':
                    self.assertEqual(sim.comp1ROI.vs1R1['fwd'].Extent, sim.comp1.vs1R1['fwd'].Extent)
                    sim.comp1ROI.Diffc1S1.Extent
                    # self.assertEqual(sim.comp1ROI.Diffc1S1.Extent, sim.comp1.Diffc1S1.Extent)

                self.assertEqual(sim.patchROI.Ex.Count, sim.patch.Ex.Count)
                self.assertEqual(sim.patchROI.Ex.Amount, sim.patch.Ex.Amount)
                if sim._solverStr != 'TetOpSplit':
                    self.assertEqual(sim.patchROI.ss1R3['fwd'].Extent, sim.patch.ss1R3['fwd'].Extent)

                self.assertEqual(sim.comp1ROI[:].S1.Count, sim.TETS(self.newGeom.comp1.tets).S1.Count)
                self.assertEqual(sim.patchROI[:].Ex.Count, sim.TRIS(self.newGeom.patch.tris).Ex.Count)

        with self.assertRaises(SimPathInvalidPath):
            sim.comp1.S1['test'].Count

    def _test_API2_SimPathSetSyntax(self, sim, usingMesh=False):
        self.init_API2_sim(sim)
    
        if usingMesh:
            self.c1DS1 *= 0.75
            self.c1DS2 *= 1.234
            self.c2DS1 *= 2.654
            self.c2DS2 *= 0.124
        else:
            self.v1 *= 2
            self.v2 *= 1.324
            self.a1 *= 3
            sim.comp1.Vol = self.v1
            sim.comp2.Vol = Params(self.v2)
            sim.patch.Area = self.a1

        self.r01f *= 1.5
        self.r01b /= 1.5
        self.r02f *= 1.3
        self.r02b = self.r02f
        self.r03f *= 3.2
        self.r03b *= 2.75
        self.r04f *= 6.421
        self.r04b = self.r04f
        self.r05f *= 1.75
        self.r05b *= 3.82
        self.cr01f *= 1.256
        self.cr01b *= 0.562
        self.initC1S1 = int(self.initC1S1 * 2.456)
        self.initP1Ex = int(self.initP1Ex * 1.658)
        self.initC1S2 = int(self.initC1S2 * 1.657)
        self.initC2S1 = int(self.initC2S1 * 1.584)
        self.initP1ExS1 = int(self.initP1ExS1 * 3.214)

        # Solver global values
        if usingMesh:
            sim.Temp = self.solvTemp * 1.245
            self.assertAlmostEqualWThresh(sim.Temp, self.solvTemp * 1.245, self.tolerance)
            sim.Temp = self.solvTemp

            with self.assertRaises(AttributeError):
                sim.temp = self.solvTemp

        # Single value setting
        if not self.useDist:
            sim.comp1.vs1R1['fwd'].K = self.r01f
            sim.comp1.vs1R1['bkw'].K = self.r01b
            sim.patch.ss1R3['fwd'].K = self.r03f
            sim.patch.ss1R3['bkw'].K = self.r03b
            sim.comp1.vs1CR1['fwd'].K = self.cr01f
            sim.comp1.vs1CR1['bkw'].K = self.cr01b
            sim.comp1.vs1R1['fwd'].Active = False
            sim.patch.ss1R5['fwd'].Active = False
            self.assertEqual(sim.comp1.vs1R1['fwd'].Active, False)
            self.assertEqual(sim.comp1.vs1R1['bkw'].Active, True)
            self.assertEqual(sim.patch.ss1R5['fwd'].Active, False)
            self.assertEqual(sim.patch.ss1R5['bkw'].Active, True)
        sim.comp1.S1.Count = self.initC1S1
        sim.patch.Ex.Count = self.initP1Ex
        if not self.useDist:
            sim.comp1.S2.Conc = self.initC1S2 / Avogad / (1e3 * sim.comp1.Vol)
            sim.comp2.S1.Amount = self.initC2S1 / Avogad
            sim.patch.ExS1.Amount = self.initP1ExS1 / Avogad
            sim.comp1.S1.Clamped = True
            sim.patch.Ex.Clamped = True
            self.assertTrue(sim.comp1.S1.Clamped)
            self.assertTrue(sim.patch.Ex.Clamped)

        if usingMesh:
            sus1, sus2 = self.newMdl.sus1, self.newMdl.sus2
            if not self.useDist:
                sim.comp1.Diffc1S1.D = self.c1DS1
                sim.comp1.Diffc1S2.D = self.c1DS2
                sim.comp2.Diffc2S1.D = self.c2DS1
                sim.comp2.Diffc2S2.D = self.c2DS2
                sim.comp1.Diffc1CC[sus1, sus1].D = 123
                sim.comp1.Diffc1CC[sus1, sus2].D = 456
                sim.comp1.Diffc1CC[sus2, sus2].D = 789
                self.assertEqual(sim.comp1.Diffc1CC[sus1, sus1].D, 123)
                self.assertEqual(sim.comp1.Diffc1CC[sus1, sus2].D, 456)
                self.assertEqual(sim.comp1.Diffc1CC[sus2, sus2].D, 789)
                self.c1DCC = CompDepDcst(lambda s: s.Count(sus1) * 1.21e-12, self.newMdl.CC)
                sim.comp1.Diffc1CC.D = self.c1DCC

                sim.comp1.Diffc1S1.Active = False
                self.assertFalse(sim.comp1.Diffc1S1.Active)
                sim.comp1.Diffc1S1.Active = True

            sim.diffb.S1.DiffusionActive = True
            self.assertTrue(sim.diffb.S1.DiffusionActive)
            if not self.useDist:
                sim.diffb.S1.Dcst = self.c1DS1 * 1.2335
                # self.assertAlmostEqualWThresh(sim.diffb.S1.Dcst, self.c1DS1 * 1.2335, self.tolerance)
                sim.diffb.S1.Dcst = self.c1DS1
            sim.diffb.S1.DiffusionActive = False

            if not self.useDist:
                sim.sdiffb.Ex.DiffusionActive = True
                self.assertTrue(sim.sdiffb.Ex.DiffusionActive)
                sim.sdiffb.Ex.Dcst = self.patchDEx * 1.2335
                # self.assertAlmostEqualWThresh(sim.sdiffb.Ex.Dcst, self.patchDEx * 1.2335, self.tolerance)
                sim.sdiffb.Ex.Dcst = self.patchDEx
                sim.sdiffb.Ex.DiffusionActive = False

            tri1 = self.newGeom.patch.tris[0]
            tet1 = [tet for tet in tri1.tetNeighbs if tet in self.newGeom.comp1.tets][0]
            tet2 = [tet for tet in tri1.tetNeighbs if tet in self.newGeom.comp2.tets][0]
            if self.useDist:
                vert1 = tri1.verts[0]
            else:
                vert1 = tri1.bars[0].verts[0]

            currVal = sim.TET(tet1).S1.Count
            sim.TET(tet1).S1.Count = 123
            self.assertEqual(sim.TET(tet1).S1.Count, 123)
            if not self.useDist:
                sim.TET(tet1).S1.Amount = 456 / Avogad
                self.assertAlmostEqualWThresh(sim.TET(tet1).S1.Amount, 456 / Avogad, self.tolerance)
                sim.TET(tet1).S1.Conc = 789 / Avogad / tet1.Vol
                self.assertAlmostEqualWThresh(sim.TET(tet1).S1.Conc, 789 / Avogad / tet1.Vol, self.tolerance)
                sim.TET(tet1).S1.Clamped = True
                self.assertTrue(sim.TET(tet1).S1.Clamped)
                sim.TET(tet1).S1.Clamped = False
            sim.TET(tet1).S1.Count = currVal

            if not self.useDist:
                sim.TET(tet1).vs1R1['fwd'].K = self.r01f * 2
                self.assertAlmostEqualWThresh(sim.TET(tet1).vs1R1['fwd'].K, self.r01f * 2, self.tolerance)
                sim.TET(tet1).vs1R1['fwd'].Active = False
                self.assertFalse(sim.TET(tet1).vs1R1['fwd'].Active)
                sim.TET(tet1).vs1R1['fwd'].Active = True
                sim.TET(tet1).vs1R1['fwd'].K = self.r01f

                sim.TET(tet1).Diffc1S1.D = self.c1DS1 * 2
                self.assertAlmostEqualWThresh(sim.TET(tet1).Diffc1S1.D, self.c1DS1 * 2, self.tolerance)
                sim.TET(tet1).Diffc1S1.D = self.c1DS1
                sim.TET(tet1).Diffc1S1.Active = False
                self.assertFalse(sim.TET(tet1).Diffc1S1.Active)
                sim.TET(tet1).Diffc1S1.Active = True

                sim.TET(tet1).V = self.membPot * 1.123
                self.assertAlmostEqualWThresh(sim.TET(tet1).V, self.membPot * 1.123, self.tolerance)
                sim.TET(tet1).V = self.membPot
                sim.TET(tet1).VClamped = True
                self.assertTrue(sim.TET(tet1).VClamped)
                sim.TET(tet1).VClamped = False

            currVal = sim.TRI(tri1).Ex.Count
            sim.TRI(tri1).Ex.Count = 124
            self.assertEqual(sim.TRI(tri1).Ex.Count, 124)
            if not self.useDist:
                sim.TRI(tri1).Ex.Amount = 456 / Avogad
                self.assertAlmostEqualWThresh(sim.TRI(tri1).Ex.Amount, 456 / Avogad, self.tolerance)
                sim.TRI(tri1).Ex.Clamped = True
                self.assertTrue(sim.TRI(tri1).Ex.Clamped)
                sim.TRI(tri1).Ex.Clamped = False
            sim.TRI(tri1).Ex.Count = currVal

            if not self.useDist:
                sim.TRI(tri1).ss1R3['fwd'].K = self.r03f * 2
                self.assertAlmostEqualWThresh(sim.TRI(tri1).ss1R3['fwd'].K, self.r03f * 2, self.tolerance)
                sim.TRI(tri1).ss1R3['fwd'].Active = False
                self.assertFalse(sim.TRI(tri1).ss1R3['fwd'].Active)
                sim.TRI(tri1).ss1R3['fwd'].Active = True
                sim.TRI(tri1).ss1R3['fwd'].K = self.r03f

                sim.TRI(tri1).DiffpatchEx.D = self.patchDEx * 2
                self.assertAlmostEqualWThresh(sim.TRI(tri1).DiffpatchEx.D, self.patchDEx * 2, self.tolerance)
                sim.TRI(tri1).DiffpatchEx.D = self.patchDEx
                # sim.TRI(tri1).DiffpatchEx.Active = False
                # self.assertFalse(sim.TRI(tri1).DiffpatchEx.Active)
                # sim.TRI(tri1).DiffpatchEx.Active = True

                sim.TRI(tri1).V = self.membPot * 1.598
                self.assertAlmostEqualWThresh(sim.TRI(tri1).V, self.membPot * 1.598, self.tolerance)
                sim.TRI(tri1).V = self.membPot
                sim.TRI(tri1).VClamped = True
                self.assertTrue(sim.TRI(tri1).VClamped)
                sim.TRI(tri1).VClamped = False

                sim.TRI(tri1).IClamp = 1.234e-7
                self.assertEqual(sim.TRI(tri1).IClamp, 1.234e-7)
                sim.TRI(tri1).IClamp = 0

                sim.TRI(tri1).ss1RVDep01['fwd'].Active = False
                self.assertFalse(sim.TRI(tri1).ss1RVDep01['fwd'].Active)
                sim.TRI(tri1).ss1RVDep01['fwd'].Active = True

                sim.TRI(tri1).Capac = self.membCap * 1.254
                # self.assertAlmostEqualWThresh(sim.TRI(tri1).Capac, self.membCap * 1.254, self.tolerance)
            
                sim.VERT(vert1).V = self.membPot * 1.2468
                self.assertAlmostEqualWThresh(sim.VERT(vert1).V, self.membPot * 1.2468, self.tolerance)
                sim.VERT(vert1).V = self.membPot
                sim.VERT(vert1).VClamped = True
                self.assertTrue(sim.VERT(vert1).VClamped)
                sim.VERT(vert1).VClamped = False
            sim.VERT(vert1).IClamp = 1.7852e-7
            # TODO Remove condition after changes from vesicle branch have been merged
            IclampVal = sim.VERT(vert1).IClamp
            if MPI.rank == 0 or not self.useDist:
                self.assertAlmostEqualWThresh(IclampVal, 1.7852e-7, self.tolerance)
            sim.VERT(vert1).IClamp = 0

            if not self.useDist:
                # ROI
                currVal = sim.comp1ROI.S1.Count
                sim.comp1ROI.S1.Count = 452
                self.assertEqual(sim.comp1ROI.S1.Count, 452)
                sim.comp1ROI.S1.Amount = 148 / Avogad
                self.assertAlmostEqualWThresh(sim.comp1ROI.S1.Amount, 148 / Avogad, self.tolerance)
                sim.comp1ROI.S1.Conc = 785 / Avogad / sim.comp1ROI.Vol
                self.assertAlmostEqualWThresh(sim.comp1ROI.S1.Conc, 785 / Avogad / sim.comp1ROI.Vol, self.tolerance)
                sim.comp1ROI.S1.Clamped = True
                # self.assertTrue(sim.comp1ROI.S1.Clamped)
                sim.comp1ROI.S1.Clamped = False
                sim.comp1ROI.S1.Count = currVal

                sim.comp1ROI.vs1R1['fwd'].K = self.r01f * 2
                # self.assertAlmostEqualWThresh(sim.comp1ROI.vs1R1['fwd'].K, self.r01f * 2, self.tolerance)
                sim.comp1ROI.vs1R1['fwd'].Active = False
                # self.assertFalse(sim.comp1ROI.vs1R1['fwd'].Active)
                sim.comp1ROI.vs1R1['fwd'].Active = True
                sim.comp1ROI.vs1R1['fwd'].K = self.r01f

                sim.patchROI.ss1R3['fwd'].K = self.r03f * 2
                # self.assertAlmostEqualWThresh(sim.patchROI.ss1R3['fwd'].K, self.r03f * 2, self.tolerance)
                sim.patchROI.ss1R3['fwd'].Active = False
                # self.assertFalse(sim.patchROI.ss1R3['fwd'].Active)
                sim.patchROI.ss1R3['fwd'].Active = True
                sim.patchROI.ss1R3['fwd'].K = self.r03f

                sim.comp1ROI.Diffc1S1.D = self.c1DS1 * 2
                # self.assertAlmostEqualWThresh(sim.comp1ROI.Diffc1S1.D, self.c1DS1 * 2, self.tolerance)
                sim.comp1ROI.Diffc1S1.D = self.c1DS1
                sim.comp1ROI.Diffc1S1.Active = False
                # self.assertFalse(sim.comp1ROI.Diffc1S1.Active)
                sim.comp1ROI.Diffc1S1.Active = True

                sim.patchROI.ss1RVDep01['fwd'].Active = False
                # self.assertFalse(sim.patchROI.ss1RVDep01['fwd'].Active)
                sim.patchROI.ss1RVDep01['fwd'].Active = True


        if not self.useDist:
            # Setting several values
            sim.comp2.vs2R2.K = self.r02f
            sim.comp2.vs2R2.K = [self.r02f] * 2
            sim.patch.ss1R4.K = self.r04f
            sim.patch.ss1R4.K = [self.r04f] * 2
            sim.patch.ss1R5.K = [self.r05f, self.r05b]
            sim.patch.ss1R4.Active = False
            sim.patch.ss1R6.Active = [True, False]
            self.assertEqual(sim.patch.ss1R4['fwd'].Active, False)
            self.assertEqual(sim.patch.ss1R4['bkw'].Active, False)
            self.assertEqual(sim.patch.ss1R6['fwd'].Active, True)
            self.assertEqual(sim.patch.ss1R6['bkw'].Active, False)
            sim.comp1.ALL(Species).Clamped = [True, False, False] if usingMesh else [True, False]
            sim.patch.ALL(Species).Clamped = [True, False, True, False]

            sim.comp1.vs1CR1['fwd'].K = CompDepRate(lambda s: s.Count(self.newMdl.sus1) * 1234, [self.newMdl.CC])
            self.assertEqual(set(sim.comp1.vs1CR1['fwd'].K), set([1234, 1234*4]))
            sim.comp1.vs1CR1['fwd'].K = self.cr01f

            # Setting several values using special methods
            sim.comp1.ALL(Species).Clamped = False
        
            if usingMesh:
                sim.membrane.Res = Params(1, self.membPot)


            sim.ALL(Compartment, Patch).ALL(Species).Clamped = False
            sim.ALL(Compartment, Patch).ALL(Reaction).Active = True
            # Check that all dependent values are correct
            self._test_API2_SimPathGetSyntax(sim, init=False, usingMesh=usingMesh)

        # Test bad syntax
        with self.assertRaises(SimPathInvalidPath):
            sim.comp4.S1.Count = 0
        with self.assertRaises(SimPathInvalidPath):
            a = sim.comp1.S1.count
        with self.assertRaises(SimPathInvalidPath):
            a = sim.comp1.vol
        with self.assertRaises(SimPathInvalidPath):
            a = sim.comp1.area

        with self.assertRaises(SimPathSolverMissingMethod):
            sim.comp1.vs1R1['fwd'].H = 10
        with self.assertRaises(SimPathSolverMissingMethod):
            sim.comp1.vs1R1['fwd'].A = 10
        with self.assertRaises(SimPathSolverMissingMethod):
            sim.comp1.vs1R1['fwd'].C = 10
        with self.assertRaises(SimPathSolverMissingMethod):
            sim.comp1.vs1R1['fwd'].Extent = 10

        # Typos to endnames
        with self.assertRaises(AttributeError):
            sim.comp1.S1.count = 10
        with self.assertRaises(AttributeError):
            sim.comp1.S1.conc = 10e-6

        with self.assertRaises(Exception):
            sim.comp1.CC.sus1.Count = 10

        with self.assertRaises(Exception):
            # patch does not contain S1
            sim.ALL().S1.Count
        with self.assertRaises(SimPathSolverMissingMethod):
            # Cannot get the count of reactions
            sim.comp1.ALL(Reaction).Count
        with self.assertRaises(SolverCallError):
            sp = sim.comp1.S1
            self.newGeom.comp1.name = 5
            sp.Count
        self.newGeom.comp1.name = 'comp1'
        with self.assertRaises(SimPathInvalidPath):
            if usingMesh:
                sim.comp1.S95.Count
            else:
                sim.comp1.ALL(Diffusion).D

        with self.assertRaises(SimPathInvalidPath):
            sim.LIST().S1.Count
        with self.assertRaises(SimPathInvalidPath):
            sim.LIST('test').S1.Count
        with self.assertRaises(SimPathInvalidPath):
            sim.MATCH('test').S1.Count

        with self.assertRaises(SimPathInvalidPath):
            sim.comp1.ALL(Species).Count = [0] * 10


class TetTestModelFramework(TestModelFramework):

    def setUp(self):
        TestModelFramework.setUp(self)

        # Model parameters

        self.initChan1Cl = 10
        self.Chan1_G = 20e-15
        self.Chan1_rev = -77e-3
        self.Chan1_P = 2.5e-15

        self.c1DS1 = 1.2e-14
        # self.c1DS1 = 1.2e-12
        self.c1DS2 = 2.5e-12
        self.c1DS3 = 0.9e-12
        self.c2DS1 = 1.9e-12
        self.c2DS2 = 2.1e-12
        self.c2DS3 = 0.8e-12
        self.patchDEx = 2.1e-13
        self.CCmult = 1e-12

        self.vrange = [-200.0e-3, 200e-3, 1e-3]

        # Geom parameters

        self.membPot = -65e-3
        self.membCap = 1e-2
        self.volRes = 1
        self.efielddt = 1e-3
        self.solvTemp = 300

        self.patchz = 2e-6

    def get_API2_Mdl(self, sas=True):
        # new model
        nmdl = super().get_API2_Mdl(sas=sas)

        r = ReactionManager()

        self.c1DCC = CompDepDcst(lambda s: s.Count(nmdl.sus2) * self.CCmult, [nmdl.CC], name='c1DCC')

        self.rvdep01f = VDepRate(lambda v: 2.0, vrange=self.vrange, name='rvdep01f')
        self.rvdep01b = VDepRate(lambda v: 1.0, vrange=self.vrange, name='rvdep01b')

        with nmdl:

            chanop1, chanop2, chancl = SubUnitState.Create()
            Chan1 = Channel.Create([chanop1, chanop2, chancl])

            nmdl.S3.valence = 1
            with nmdl.ssys, Chan1[...]:
                chancl.s <r['ss1RVDep01']> chanop1.s
                r['ss1RVDep01'].K = self.rvdep01f, self.rvdep01b

                Chan1_Ohm_I = OhmicCurr.Create(Chan1[chanop1 | chanop2], self.Chan1_G, self.Chan1_rev)
                Chan1_GHK_I = GHKCurr.Create(Chan1[chanop1], nmdl.S3, self.Chan1_P, computeflux=True)

            with nmdl.vsys1:
                Diffc1S1, Diffc1S2, Diffc1S3, Diffc1CC = Diffusion.Create(
                    Params(nmdl.S1, self.c1DS1),
                    Params(nmdl.S2, self.c1DS2),
                    Params(nmdl.S3, self.c1DS3),
                    Params(nmdl.CC, self.c1DCC),
                )
            with nmdl.vsys2:
                Diffc2S1, Diffc2S2, Diffc2S3 = Diffusion.Create(
                    Params(nmdl.S1, self.c2DS1),
                    Params(nmdl.S2, self.c2DS2),
                    Params(nmdl.S3, self.c2DS3),
                )

            if not self.useDist:
                with nmdl.ssys:
                    DiffpatchEx = Diffusion.Create(nmdl.Ex, self.patchDEx)

        return nmdl

    def get_API1_Mdl(self):
        # old model
        omdl = super().get_API1_Mdl()

        omdl.getSpec('S3').setValence(1)
        Chan1 = smodel.Chan('Chan1', omdl)
        chanop1 = smodel.ChanState('chanop1', omdl, Chan1)
        chanop2 = smodel.ChanState('chanop2', omdl, Chan1)
        chancl = smodel.ChanState('chancl', omdl, Chan1)
        ss1RVDep01_fwd = smodel.VDepSReac('ss1RVDep01_fwd', omdl.getSurfsys('ssys'), slhs=[chancl],
                srhs=[chanop1], k=self.rvdep01f._func, vrange=self.rvdep01f.vrange)
        ss1RVDep01_bkw = smodel.VDepSReac('ss1RVDep01_bkw', omdl.getSurfsys('ssys'), srhs=[chancl],
                slhs=[chanop1], k=self.rvdep01b._func, vrange=self.rvdep01b.vrange)

        Chan1_Ohm_I_1 = smodel.OhmicCurr('Chan1_Ohm_I_1', omdl.getSurfsys('ssys'), chanstate=chanop1, g=self.Chan1_G, erev=self.Chan1_rev)
        Chan1_Ohm_I_2 = smodel.OhmicCurr('Chan1_Ohm_I_2', omdl.getSurfsys('ssys'), chanstate=chanop2, g=self.Chan1_G, erev=self.Chan1_rev)
        Chan1_GHK_I = smodel.GHKcurr('Chan1_GHK_I', omdl.getSurfsys('ssys'), chanop1, omdl.getSpec('S3'), computeflux = True)
        Chan1_GHK_I.setP(self.Chan1_P)

        Diffc1S1 = smodel.Diff('Diffc1S1', omdl.getVolsys('vsys1'), omdl.getSpec('S1'))
        Diffc1S1.setDcst(self.c1DS1)
        Diffc1S2 = smodel.Diff('Diffc1S2', omdl.getVolsys('vsys1'), omdl.getSpec('S2'))
        Diffc1S2.setDcst(self.c1DS2)
        Diffc1S3 = smodel.Diff('Diffc1S3', omdl.getVolsys('vsys1'), omdl.getSpec('S3'))
        Diffc1S3.setDcst(self.c1DS3)
        Diffc1CC11 = smodel.Diff('Diffc1CC11', omdl.getVolsys('vsys1'), omdl.getSpec('CCsus11'))
        Diffc1CC11.setDcst(0)
        Diffc1CC12 = smodel.Diff('Diffc1CC12', omdl.getVolsys('vsys1'), omdl.getSpec('CCsus12'))
        Diffc1CC12.setDcst(self.CCmult)
        Diffc1CC22 = smodel.Diff('Diffc1CC22', omdl.getVolsys('vsys1'), omdl.getSpec('CCsus22'))
        Diffc1CC22.setDcst(2 * self.CCmult)

        Diffc2S1 = smodel.Diff('Diffc2S1', omdl.getVolsys('vsys2'), omdl.getSpec('S1'))
        Diffc2S1.setDcst(self.c2DS1)
        Diffc2S2 = smodel.Diff('Diffc2S2', omdl.getVolsys('vsys2'), omdl.getSpec('S2'))
        Diffc2S2.setDcst(self.c2DS2)
        Diffc2S3 = smodel.Diff('Diffc2S3', omdl.getVolsys('vsys2'), omdl.getSpec('S3'))
        Diffc2S3.setDcst(self.c2DS3)

        if not self.useDist:
            DiffpatchEx = smodel.Diff('DiffpatchEx', omdl.getSurfsys('ssys'), omdl.getSpec('Ex'))
            DiffpatchEx.setDcst(self.patchDEx)

        return omdl

    def get_API2_MeshOnly(self):
        return TetMesh.Load(os.path.join(FILEDIR, '../geom_test/meshes', 'cyl_len10_diam1'), name='nmesh')

    def get_API2_Geom(self, nmdl):

        nmesh = self.get_API2_MeshOnly()

        with nmesh:
            compKwArgs = dict(conductivity=1 / self.volRes) if self.useDist else {}

            center1 = nmesh.tets[0, 0, nmesh.bbox.max.z / 2]
            self.c1tets = TetList([tet for tet in nmesh.tets if tet.center.z > self.patchz and tet.idx != center1.idx])
            c2tets = nmesh.tets - self.c1tets - TetList([center1])
            comp1 = Compartment.Create(self.c1tets, nmdl.vsys1, **compKwArgs)
            comp1bis = Compartment.Create([center1], nmdl.vsys1, **compKwArgs)
            comp2 = Compartment.Create(c2tets, nmdl.vsys2, **compKwArgs)
            ptris = self.c1tets.surface & c2tets.surface
            self.ptris = TriList(tri for tri in ptris if tri.center.x > 0)
            self.ptris2 = ptris - self.ptris
            patch = Patch.Create(self.ptris, comp1, comp2, nmdl.ssys)
            patch2 = Patch.Create(self.ptris2, comp1, comp2, nmdl.ssys)

            if not self.useDist:
                comp1ROI = ROI.Create(comp1.tets)
                patchROI = ROI.Create(patch.tris)
                vertROI = ROI.Create(patch.tris.verts)

                membrane = Membrane.Create([patch, patch2], opt_method = 1)

                sdiffb = SDiffBoundary.Create(patch.tris.edges & patch2.tris.edges, patch, patch2)
            else:
                membrane = Membrane.Create([patch], capacitance=self.membCap)
                membrane2 = Membrane.Create([patch2], capacitance=self.membCap)

            diffb = DiffBoundary.Create(center1.faces)

        self.v1 = comp1.Vol
        self.v2 = comp2.Vol
        self.a1 = patch.Area

        return nmesh

    def get_API1_MeshOnly(self):
        return smeshio.loadMesh(os.path.join(FILEDIR, '../geom_test/meshes', 'cyl_len10_diam1'))[0]

    def get_API1_Geom(self, omdl):

        omesh = self.get_API1_MeshOnly()

        n = omesh.countTets()
        center1 = omesh.findTetByPoint([0, 0, omesh.getBoundMax()[2] / 2])
        self.oc1tets = [i for i in range(n) if omesh.getTetBarycenter(i)[2] > self.patchz and i != center1]
        oc1btets = [center1]
        oc2tets = [i for i in range(n) if omesh.getTetBarycenter(i)[2] <= self.patchz and i != center1]

        trioverlap = smeshctrl.findOverlapTris(omesh, self.oc1tets, oc2tets)
        self.optris = [i for i in trioverlap if omesh.getTriBarycenter(i)[0] > 0]
        self.optris2 = [i for i in trioverlap if omesh.getTriBarycenter(i)[0] <= 0]

        comp1 = sgeom.TmComp('comp1', omesh, self.oc1tets)
        comp1.addVolsys('vsys1')
        comp1bis = sgeom.TmComp('comp1bis', omesh, oc1btets)
        comp1.addVolsys('vsys1')
        comp2 = sgeom.TmComp('comp2', omesh, oc2tets)
        comp2.addVolsys('vsys2')
        patch = sgeom.TmPatch('patch', omesh, self.optris, comp1, comp2)
        patch.addSurfsys('ssys')
        patch2 = sgeom.TmPatch('patch2', omesh, self.optris2, comp1, comp2)
        patch2.addSurfsys('ssys')

        if not self.useDist:
            omesh.addROI('comp1ROI', sgeom.ELEM_TET, self.oc1tets)
            omesh.addROI('patchROI', sgeom.ELEM_TRI, self.optris)
            verts = []
            for i in self.optris:
                verts += omesh.getTri(i)
            omesh.addROI('vertROI', sgeom.ELEM_VERTEX, verts)

            bars1 = set()
            for i in self.optris:
                bars1 |= set(omesh.getTriBars(i))
            bars2 = set()
            for i in self.optris2:
                bars2 |= set(omesh.getTriBars(i))
            bars = list(bars1 & bars2)
            sdiffb = sgeom.SDiffBoundary('sdiffb', omesh, bars, [patch, patch2])

        membrane = sgeom.Memb('membrane', omesh, [patch, patch2], opt_method = 1)

        diffb = sgeom.DiffBoundary('diffb', omesh, omesh.getTetTriNeighb(center1))

        self.v1 = comp1.getVol()
        self.v2 = comp2.getVol()
        self.a1 = patch.getArea()

        return omesh

    def init_API2_sim(self, sim):
        super().init_API2_sim(sim)
        sim.comp1.S3.Count = self.initC1S3
        sim.comp2.S3.Count = self.initC2S3
        sim.membrane.Potential = self.membPot
        if not self.useDist:
            sim.membrane.Capac = self.membCap
            sim.membrane.VolRes = self.volRes
        else:
            sim.membrane2.Potential = self.membPot
        sim.Temp = self.solvTemp

        sim.patch.Chan1[sim.model.chancl].Count = self.initChan1Cl

    def init_API1_sim(self, sim):
        super().init_API1_sim(sim)
        sim.setCompCount('comp1', 'S3', self.initC1S3)
        sim.setCompCount('comp2', 'S3', self.initC2S3)
        sim.setMembPotential('membrane', self.membPot)
        sim.setMembCapac('membrane', self.membCap)
        sim.setMembVolRes('membrane', self.volRes)
        sim.setTemp(self.solvTemp)

        sim.setPatchCount('patch', 'chancl', self.initChan1Cl)

    def _test_API2_SimPathSetSyntax(self, sim, usingMesh=True):
        super()._test_API2_SimPathSetSyntax(sim, usingMesh=True)

    def _test_API2_SimPathGetSyntax(self, sim, init=True, usingMesh=True, extents=None):
        super()._test_API2_SimPathGetSyntax(sim, init=init, usingMesh=True, extents=extents)


class DistTetTestModelFramework(TetTestModelFramework):

    def setUp(self):
        TetTestModelFramework.setUp(self)
        self.useDist = True
        self.endTime /= 10

    def get_API2_MeshOnly(self):
        return DistMesh(os.path.join(FILEDIR, '../geom_test/meshes/cyl_len10_diam1.msh'), scale=1e-6, name='nmesh')

    def get_API1_MeshOnly(self):
        return smeshio.importGmsh(os.path.join(FILEDIR, '../geom_test/meshes/cyl_len10_diam1.msh'), 1e-6)[0]

