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

""" Unit tests for data saving."""

import os
import tempfile
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

from . import test_model


class WmCheckpoints(test_model.TestModelFramework):
    """Test checkpoint and restore methods for Wmdirect solver"""
    def setUp(self, callParent=True):
        if callParent:
            super().setUp()

        self.newMdl = self.get_API2_Mdl()
        self.newGeom = self.get_API2_Geom(self.newMdl)

        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)

        self.cpTime = self.endTime / 10
        self.nbRuns = 2

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Wmdirect', nmdl, ngeom, nrng)
        return nsim

    def _setInitValsToCurrentOnes(self):
        sim = self.newSim

        sus1, sus2 = self.sus1, self.sus2

        self.initC1S1 = sim.comp1.S1.Count
        self.initC1S2 = sim.comp1.S2.Count
        self.initC2S1 = sim.comp2.S1.Count
        self.initC2S2 = sim.comp2.S2.Count
        self.initP1Ex = sim.patch.Ex.Count
        self.initP1ExS1 = sim.patch.ExS1.Count
        self.initP1ExS2 = sim.patch.ExS2.Count
        self.initP1ExS1S2 = sim.patch.ExS1S2.Count

        self.initC1CCsus11 = sim.comp1.CC[sus1, sus1].Count
        self.initC1CCsus12 = sim.comp1.CC[sus1, sus2].Count
        self.initC1CCsus22 = sim.comp1.CC[sus2, sus2].Count

        # TODO This should be set with the current extent value but Wmdirect does not save them
        # in the current STEPS version.
        self.extents = [0]*13


    def testBasicCheckpoints(self):
        cpPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}basic')
        cpEmptyPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}empty')

        self.newSim.newRun()
        self.init_API2_sim(self.newSim)

        self.newSim.run(self.endTime / 2)

        with self.assertRaises(Exception):
            self.checkpoint(42)

        self.newSim.checkpoint(cpPath)
        self._setInitValsToCurrentOnes()

        self.newSim.newRun()

        with self.assertRaises(Exception):
            self.newSim.restore(cpEmptyPath)

        self.newSim.restore(cpPath)

        self.assertAlmostEqual(self.newSim.Time, self.endTime / 2)
        self._test_API2_SimPathGetSyntax(self.newSim, init=False, extents=self.extents)

        if MPI._shouldWrite:
            os.remove(cpPath)

    def testAutoCheckpoints(self):
        tmpDir = tempfile.gettempdir()
        prefix = f'{self.__class__.__name__}auto'
        pathPrefix = os.path.join(tmpDir, prefix)

        with self.assertRaises(TypeError):
            self.newSim.autoCheckpoint(self.cpTime, 42)
        
        with self.assertRaises(ValueError):
            self.newSim.autoCheckpoint(-self.cpTime, pathPrefix)

        with self.assertRaises(ValueError):
            self.newSim.autoCheckpoint(0, pathPrefix)

        self.newSim.autoCheckpoint(self.cpTime, pathPrefix)

        for r in range(self.nbRuns):
            self.newSim.newRun()
            self.newSim.run(self.endTime)

        cpfiles = []
        for name in os.listdir(tmpDir):
            if name.startswith(prefix):
                run, time, solver, cp = name[(len(prefix)+1):].split('_')
                run = int(run)
                time = float(time)
                self.assertEqual(cp, 'cp')
                cpfiles.append((run, time, name))
        cpfiles.sort()

        self.assertEqual(len(cpfiles), self.nbRuns * (self.endTime / self.cpTime + 1))
        for r in range(self.nbRuns):
            self.assertEqual(sum(1 for run, time, _ in cpfiles if run == r), self.endTime / self.cpTime + 1)

        for run, time, name in cpfiles:
            cpPath = os.path.join(tmpDir, name)

            self.newSim.newRun()

            self.newSim.restore(cpPath)

            self.assertAlmostEqual(self.newSim.Time, time)

            if MPI._shouldWrite:
                os.remove(cpPath)

    def testOnlyLast(self):
        tmpDir = tempfile.gettempdir()
        prefix = f'{self.__class__.__name__}onlyLast'
        pathPrefix = os.path.join(tmpDir, prefix)

        self.newSim.autoCheckpoint(self.cpTime, pathPrefix, onlyLast=True)

        for r in range(self.nbRuns):
            self.newSim.newRun()
            self.newSim.run(self.endTime)

        cpfiles = []
        for name in os.listdir(tmpDir):
            if name.startswith(prefix):
                run, time, solver, cp = name[(len(prefix)+1):].split('_')
                run = int(run)
                time = float(time)
                self.assertEqual(cp, 'cp')
                cpfiles.append((run, time, name))
        cpfiles.sort()
        
        self.assertEqual(len(cpfiles), 1)
        self.assertEqual(cpfiles[0][0], self.nbRuns - 1)

        cpPath = os.path.join(tmpDir, cpfiles[0][2])

        self.newSim.newRun()

        self.newSim.restore(cpPath)

        self.assertAlmostEqual(self.newSim.Time, cpfiles[0][1])

        if MPI._shouldWrite:
            os.remove(cpPath)
        
    def testStopAutoCheckpointing(self):
        tmpDir = tempfile.gettempdir()
        prefix = f'{self.__class__.__name__}stopAutoChkpt'
        pathPrefix = os.path.join(tmpDir, prefix)

        for r in range(self.nbRuns):
            self.newSim.newRun()
            self.newSim.autoCheckpoint(self.cpTime, pathPrefix)
            self.newSim.run(self.endTime / 2)
            self.newSim.autoCheckpoint(None)
            self.newSim.run(self.endTime)

        cpfiles = []
        for name in os.listdir(tmpDir):
            if name.startswith(prefix):
                run, time, solver, cp = name[(len(prefix)+1):].split('_')
                run = int(run)
                time = float(time)
                self.assertEqual(cp, 'cp')
                cpfiles.append((run, time, name))
        cpfiles.sort()

        self.assertEqual(len(cpfiles), self.nbRuns * (self.endTime / 2 / self.cpTime + 1))
        for r in range(self.nbRuns):
            self.assertEqual(sum(1 for run, time, _ in cpfiles if run == r), self.endTime / 2 / self.cpTime + 1)
            self.assertTrue(all(time <= self.endTime / 2 for run, time, _ in cpfiles if run == r))

        for run, time, name in cpfiles:
            cpPath = os.path.join(tmpDir, name)

            self.newSim.newRun()

            self.newSim.restore(cpPath)

            self.assertAlmostEqual(self.newSim.Time, time)

            if MPI._shouldWrite:
                os.remove(cpPath)

    def testContinueRestoredRun(self):
        tmpDir = tempfile.gettempdir()
        prefix = f'{self.__class__.__name__}ContinueRestoredRun'
        pathPrefix = os.path.join(tmpDir, prefix)

        rs = ResultSelector(self.newSim)
        saver = rs.comp1.S2.Count

        self.newSim.toSave(saver, dt=self.deltaT)

        self.newSim.newRun()

        self.init_API2_sim(self.newSim)
        self.newSim.run(self.endTime / 2)
        self.newSim.checkpoint(pathPrefix)
        self.newSim.run(self.endTime)

        # Continue run in same simulation
        self.newSim.newRun()
        self.newSim.restore(pathPrefix)
        self.newSim.run(self.endTime)

        if MPI._shouldWrite:
            self.assertEqual(len(saver.time), 2)
            self.assertEqual(len(saver.data), 2)

            self.assertLess(abs(len(saver.time[1]) - (self.endTime / 2 / self.deltaT + 1)), 1.5)
            self.assertTrue(all(t >= self.endTime / 2 for t in saver.time[1]))

            n = len(saver.data[0]) - len(saver.data[1])
            firstRunData = saver.data[0:1, n:, 0:1]
            self.assertSameData(firstRunData, saver.data[1:2, :, 0:1], [self.countRefVal])

        # continue run in new simulation
        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)

        rs = ResultSelector(self.newSim)
        saver = rs.comp1.S2.Count

        self.newSim.toSave(saver, dt=self.deltaT)
        self.newSim.newRun()
        self.newSim.restore(pathPrefix)
        self.newSim.run(self.endTime)

        if MPI._shouldWrite:
            self.assertEqual(len(saver.time), 1)
            self.assertEqual(len(saver.data), 1)

            self.assertLess(abs(len(saver.time[0]) - (self.endTime / 2 / self.deltaT + 1)), 1.5)
            self.assertTrue(all(t >= self.endTime / 2 for t in saver.time[0]))

            self.assertSameData(firstRunData, saver.data[0:1, :, 0:1], [self.countRefVal])

        if MPI._shouldWrite:
            os.remove(pathPrefix)


class TetCheckpoints(test_model.TetTestModelFramework, WmCheckpoints):
    """Test checkpoint and restore methods for Tetexact solver"""

    def setUp(self):
        test_model.TetTestModelFramework.setUp(self)
        WmCheckpoints.setUp(self)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Tetexact', nmdl, ngeom, nrng, True)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _setInitValsToCurrentOnes(self):
        WmCheckpoints._setInitValsToCurrentOnes(self)

        sim = self.newSim

        tri1 = self.newGeom.patch.tris[0]
        tet1 = [tet for tet in tri1.tetNeighbs if tet in self.newGeom.comp1.tets][0]
        self.membPot = sim.TET(tet1).V

        self.extents = [
            sim.comp1.vs1R1['fwd'].Extent,
            sim.comp1.vs1R1['bkw'].Extent,
            sim.comp2.vs2R2['fwd'].Extent,
            sim.comp2.vs2R2['bkw'].Extent,
            sim.patch.ss1R3['fwd'].Extent,
            sim.patch.ss1R3['bkw'].Extent,
            sim.patch.ss1R4['fwd'].Extent,
            sim.patch.ss1R4['bkw'].Extent,
            sim.patch.ss1R5['fwd'].Extent,
            sim.patch.ss1R5['bkw'].Extent,
            sim.patch.ss1R6['fwd'].Extent,
            sim.patch.ss1R6['bkw'].Extent,
            sim.patch.ss1R7.Extent,
        ]


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(WmCheckpoints, "test"))
    all_tests.append(unittest.makeSuite(TetCheckpoints, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

