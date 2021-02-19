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

import numpy as np
import os
import tempfile
import threading
import time
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

import steps.API_1.model as smodel
import steps.API_1.geom as sgeom
import steps.API_1.rng as srng
import steps.API_1.solver as ssolver
import steps.API_1.utilities.meshio as smeshio

from . import test_model

class SimDataSaving(test_model.TestModelFramework):
    """Test data access, setting, and saving."""
    def setUp(self, callParent=True):
        if callParent:
            super().setUp()

        self.newMdl = self.get_API2_Mdl()
        self.oldMdl = self.get_API1_Mdl()
        self.newGeom = self.get_API2_Geom(self.newMdl)
        self.oldGeom = self.get_API1_Geom(self.oldMdl)

        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)
        self.oldSim = self._get_API1_Sim(self.oldMdl, self.oldGeom)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Wmdirect', nmdl, ngeom, nrng)
        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        orng=srng.create('mt19937',512)
        orng.initialize(self.seed)
        osim = ssolver.Wmdirect(omdl, ogeom, orng)
        return osim

    def testSimPathGetSyntax(self, init=True, usingMesh=False):
        self._test_API2_SimPathGetSyntax(self.newSim, init=init, usingMesh=usingMesh)

    def testSimPathSetSyntax(self, usingMesh=False):
        self._test_API2_SimPathSetSyntax(self.newSim, usingMesh=usingMesh)

    def testResultSelectorSyntax(self):
        with self.assertRaises(Exception):
            ResultSelector(42)

        ns = self._get_API2_Sim(self.newMdl, self.newGeom)

        rs = ResultSelector(self.newSim)
        rs2 = ResultSelector(ns)

        # Valid
        rs.comp1.S1.Count
        rs.comp1.S1.Conc

        rs2.comp1.S1.Count

        2 * rs.comp1.S1.Count
        rs.comp1.S1.Count * 2
        2 + rs.comp1.S1.Count
        rs.comp1.S1.Count + 2
        2 / rs.comp1.S1.Count
        rs.comp1.S1.Count / 2
        2 - rs.comp1.S1.Count
        rs.comp1.S1.Count - 2
        rs.comp1.S1.Count ** 2

        rs.comp1.LIST('S1', 'S2').Count + rs.comp2.S1.Count
        rs.comp2.S1.Count + rs.comp1.LIST('S1', 'S2').Count
        rs.comp1.LIST('S1', 'S2').Count - rs.comp2.S1.Count
        rs.comp1.LIST('S1', 'S2').Count / rs.comp2.S1.Count
        rs.comp1.LIST('S1', 'S2').Count * rs.comp2.S1.Count
        rs.comp2.S1.Count * rs.comp1.LIST('S1', 'S2').Count 

        # TODO write code to test this
        rs.SUM(rs.comp1.ALL(Species).Count)

        #Invalid
        with self.assertRaises(Exception):
            rs.test

        with self.assertRaises(Exception):
            rs.comp1.S1.count

        with self.assertRaises(Exception):
            a = rs.comp1.S1
            self.newSim.toSave(a, dt=self.deltaT)

        with self.assertRaises(Exception):
            rs.comp1.S1 << rs.comp2.S2.Count
        with self.assertRaises(Exception):
            rs << rs.comp2.S2.Count

        with self.assertRaises(Exception):
            self.newSim.toSave(42)

        with self.assertRaises(TypeError):
            rs.comp1.S1.Count << 'test'
        with self.assertRaises(TypeError):
            rs.comp1.S1.Count + 'test'
        with self.assertRaises(Exception):
            rs.comp1.S1.Count << rs2.comp1.S1.Count
        with self.assertRaises(Exception):
            rs.comp1.LIST('S1', 'S2', 'CC').Count + rs.comp2.LIST('S1', 'S2').Count
        with self.assertRaises(Exception):
            rs.comp1.LIST('S1', 'S2', 'CC').Count - rs.comp2.LIST('S1', 'S2').Count
        with self.assertRaises(Exception):
            rs.comp1.LIST('S1', 'S2', 'CC').Count * rs.comp2.LIST('S1', 'S2').Count
        with self.assertRaises(Exception):
            rs.comp1.LIST('S1', 'S2', 'CC').Count / rs.comp2.LIST('S1', 'S2').Count

        with self.assertRaises(Exception):
            a = rs.comp1.S1.Count
            a.save()

    def simulateAndCompare(self, oldSimToSave, newSimToSave, refVal, dbPath=None, dbGroupArgs={}, cfreq=100):
        # Old interface
        allOldDat = []
        for rid in range(self.nbRuns):
            self.oldSim.reset()
            self.init_API1_sim(self.oldSim)

            timepoints, oldSaveFunc = oldSimToSave
            oldDat = []
            for t in timepoints:
                self.oldSim.run(t)
                oldDat.append(oldSaveFunc(self.oldSim))

            allOldDat.append(oldDat)
        allOldDat = np.array(allOldDat)

        # Delete the solver from API_1 to avoid issues with MPI
        self.oldSim = None

        # New interface
        getAndAddSaverFunc, getDataFunc, runFunc, expLabels = newSimToSave
        rs = ResultSelector(self.newSim)
        savers = getAndAddSaverFunc(rs)

        if dbPath is None:
            for rid in range(self.nbRuns):
                self.newSim.newRun()
                self.init_API2_sim(self.newSim)
                runFunc(self.newSim, savers)
            if MPI._shouldWrite:
                newDat = getDataFunc(savers)
                labelList = [s.labels for s in savers]
        else:
            # Test database saving
            with SQLiteDBHandler(dbPath, commitFreq=cfreq) as dbh:
                self.newSim.toDB(dbh, dbPath, **dbGroupArgs)
                for rid in range(self.nbRuns):
                    self.newSim.newRun()
                    self.init_API2_sim(self.newSim)
                    runFunc(self.newSim, savers)

                if MPI._shouldWrite:
                    newDat = getDataFunc(savers)
                    labelList = [s.labels for s in savers]

            # Data checking is only done in rank 0
            if MPI._shouldWrite:
                with SQLiteDBHandler(dbPath) as dbh:
                    for name, val in dbGroupArgs.items():
                        self.assertEqual(getattr(dbh[dbPath], name), val)
                    loadedSavers = dbh[dbPath].results
                    if MPI._shouldWrite:
                        newDat2 = getDataFunc(loadedSavers)
                        labelList2 = [s.labels for s in loadedSavers]

            # Check that the loaded data is identical to the one extracted during simulation
            if MPI._shouldWrite:
                self.assertSameData(newDat, newDat2, refVal)
                for llst, explbl in zip(labelList, labelList2):
                    self.assertEqual(len(llst), len(explbl))
                    for rslbl, elbl in zip(llst, explbl):
                        self.assertEqual(rslbl, elbl)

        if MPI._shouldWrite:
            # Check that the new interface gives the same results as the old one
            self.assertSameData(allOldDat, newDat, refVal)

            # Also check labels
            for llst, explbl in zip(labelList, expLabels):
                self.assertEqual(len(llst), len(explbl))
                for rslbl, elbl in zip(llst, explbl):
                    if elbl is not None:
                        self.assertEqual(rslbl, elbl)

    def _testDeltaTSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompCount('comp1', 'S1'))
            res.append(sim.getCompCount('comp1', 'S2'))
            res.append(sim.getPatchCount('patch', 'ExS1S2'))
            res.append(sim.getCompCount('comp1', 'CCsus12'))
            return res

        def newSave(rs):
            saver = \
                rs.comp1.S1.Count << \
                rs.comp1.S2.Count << \
                rs.patch.ExS1S2.Count << \
                rs.comp1.CC[self.newMdl.sus1, self.newMdl.sus2].Count
            self.newSim.toSave(saver, dt=self.deltaT)
            return [saver]

        def newGetData(saver):
            saver = saver[0]
            saver.time[0]
            np.array(saver.time)
            saver.data[0,:]
            np.array(saver.data)
            return saver.data[:]

        def newRun(sim, rs):
            sim.run(self.endTime)

        explbls = [['comp1.S1.Count', 'comp1.S2.Count', 'patch.ExS1S2.Count', 'comp1.CC_sus1_sus2.Count']]

        return (
            (timepoints, oldSave),
            (newSave, newGetData, newRun, explbls),
            [self.countRefVal]*3,
        )

    def testDeltaTSaving(self):
        self.simulateAndCompare(*self._testDeltaTSavingParams())

    def testDeltaTSavingDB(self):
        dbPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}testDeltaTSaving.db')
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testDeltaTSavingParams(), dbPath, dbdct)

    def _testTimepointSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompCount('comp1', 'S1'))
            res.append(sim.getCompCount('comp1', 'S2'))
            res.append(sim.getPatchCount('patch', 'ExS1S2'))
            return res

        def newSave(rs):
            saver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
            self.newSim.toSave(saver, timePoints=timepoints)
            return [saver]

        def newRun(sim, rs):
            sim.run(self.endTime)

        explbls = [['comp1.S1.Count', 'comp1.S2.Count', 'patch.ExS1S2.Count']]

        return (
            (timepoints, oldSave),
            (newSave, lambda s: s[0].data[:], newRun, explbls),
            [self.countRefVal]*3
        )

    def testTimepointSaving(self):
        self.simulateAndCompare(*self._testTimepointSavingParams())

    def testTimepointSavingDB(self):
        dbPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}testTimepointSaving.db')
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testTimepointSavingParams(), dbPath, dbdct)

    def _testUnspecifiedSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompCount('comp1', 'S1'))
            res.append(sim.getCompCount('comp1', 'S2'))
            res.append(sim.getPatchCount('patch', 'ExS1S2'))
            return res

        def newSave(rs):
            saver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
            self.newSim.toSave(saver)
            return [saver]

        def newRun(sim, savers):
            for t in timepoints:
                sim.run(t)
                for rs in savers:
                    rs.save()

        explbls = [['comp1.S1.Count', 'comp1.S2.Count', 'patch.ExS1S2.Count']]

        return (
            (timepoints, oldSave),
            (newSave, lambda s: s[0].data[:], newRun, explbls),
            [self.countRefVal]*3
        )

    def testUnspecifiedSaving(self):
        self.simulateAndCompare(*self._testUnspecifiedSavingParams())

    def testUnspecifiedSavingDB(self):
        dbPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}testUnspecifiedSaving.db')
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testUnspecifiedSavingParams(), dbPath, dbdct)

    def _testFileSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        path = os.path.join(tempfile.gettempdir(), 'NewData.dat')
        def oldSave(sim):
            res = []
            res.append(sim.getCompCount('comp1', 'S1'))
            res.append(sim.getCompCount('comp1', 'S2'))
            res.append(sim.getPatchCount('patch', 'ExS1S2'))
            res.append(sim.getCompCount('comp1', 'CCsus12') + 2*sim.getCompCount('comp1', 'CCsus12'))
            tmp = sim.getCompCount('comp1', 'S1') + sim.getCompCount('comp1', 'S2') + sim.getCompCount('comp2', 'S1') + sim.getCompCount('comp2', 'S2')
            tmp /= (sim.getPatchCount('patch', 'ExS1S2') + sim.getPatchCount('patch', 'Ex'))
            res.append(tmp)
            return res

        def newSave(rs):
            saver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
            saver <<= rs.comp1.CC.sus2.Count
            saver <<= rs.SUM(rs.LIST('comp1', 'comp2').LIST('S1', 'S2').Count) / rs.SUM(rs.patch.ExS1S2.Count << rs.patch.Ex.Count)
            self.newSim.toSave(saver, dt=self.deltaT)
            saver.toFile(path)
            return [saver]

        def newGetData(savers):
            saver = savers[0]
            saver.time[0]
            np.array(saver.time)
            saver.data[0,:]
            np.array(saver.data)
            rs = ResultSelector.FromFile(path)
            rs.data[0,:]
            np.array(rs.data)
            self.assertSequenceEqual(
                list(rs.data[self.nbRuns-1, 0]),
                list(rs.data[-1, 0])
            )
            self.assertEqual(self.nbRuns, len(rs.data))
            with self.assertRaises(IndexError):
                rs.data[self.nbRuns]
            with self.assertRaises(IndexError):
                rs.data[self.nbRuns:]
            return rs.data[:]

        def newRun(sim, savers):
            sim.run(self.endTime / 2)
            if MPI._shouldWrite:
                # Accessing data in the middle of a run
                for rs in savers:
                    rs.time[0]
                    rs.data[0]
            sim.run(self.endTime)

        explbls = [['comp1.S1.Count', 'comp1.S2.Count', 'patch.ExS1S2.Count', 'comp1.CC.sus2.Count', '(SUM(LIST(comp1, comp2).LIST(S1, S2).Count) / SUM(patch.ExS1S2.Count, patch.Ex.Count))']]

        return (
            (timepoints, oldSave),
            (newSave, newGetData, newRun, explbls),
            [self.countRefVal]*4
        )

    def testFileSaving(self):
        self.simulateAndCompare(*self._testFileSavingParams())

    def testFileSavingDB(self):
        dbPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}testFileSaving.db')
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testFileSavingParams(), dbPath, dbdct)

    def testPartialDataReading(self):
        path = os.path.join(tempfile.gettempdir(), 'NewData.dat')

        # Data reading function
        def readData(path):
            rs = ResultSelector.FromFile(path)
            while True:
                nbIndexErr = 0
                try:
                    rs.data[self.nbRuns-1, 0]
                except IndexError:
                    nbIndexErr += 1
                try:
                    rs.time[self.nbRuns-1, 0]
                except IndexError:
                    nbIndexErr += 1

                n = len(rs.data) - 1
                if n < self.nbRuns - 1:
                    self.assertEqual(nbIndexErr, 2)
                    try:
                        rs.data[n, -1]
                        rs.time[n, -1]
                    except IndexError:
                        pass
                else:
                    try:
                        rs.data[n, -1]
                        return None
                    except:
                        pass

        rs = ResultSelector(self.newSim)
        saver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
        saver.toFile(path)
        self.newSim.toSave(saver, dt=self.deltaT)

        if MPI._shouldWrite:
            # Reading thread
            thr = threading.Thread(target=readData, args=(path,))
            thr.start()

        # Simulation
        for rid in range(self.nbRuns):
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.endTime // 2)
            time.sleep(0.25)
            if MPI._shouldWrite:
                saver.time[0]
            time.sleep(0.25)
            self.newSim.run(self.endTime)
            time.sleep(0.25)
            if MPI._shouldWrite:
                saver.time[0]
            time.sleep(0.25)
        time.sleep(0.25)
        if MPI._shouldWrite:
            saver.time[0]

        if MPI._shouldWrite:
            thr.join(1)
            self.assertFalse(thr.is_alive())

    def _testMultiSaversParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompCount('comp1', 'S1'))
            res.append(sim.getCompCount('comp1', 'S2'))
            res.append(sim.getPatchCount('patch', 'ExS1S2'))
            res.append(sim.getCompCount('comp2', 'S1'))
            res.append(sim.getCompCount('comp2', 'S2'))
            return res

        def newSave(rs):
            saver1 = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
            saver2 = rs.comp2.S1.Count << rs.comp2.S2.Count
            self.newSim.toSave(saver1, saver2, dt=self.deltaT)
            return [saver1, saver2]

        def newGetData(savers):
            allRes = []
            for r in range(self.nbRuns):
                run = []
                for t, row in enumerate(savers[0].data[r]):
                    res = []
                    for rs in savers:
                        for val in rs.data[r, t]:
                            res.append(val)
                    run.append(res)
                allRes.append(run)
            return np.array(allRes)

        def newRun(sim, rs):
            sim.run(self.endTime)

        explbls = [['comp1.S1.Count', 'comp1.S2.Count', 'patch.ExS1S2.Count'], ['comp2.S1.Count', 'comp2.S2.Count']]

        return (
            (timepoints, oldSave),
            (newSave, newGetData, newRun, explbls),
            [self.countRefVal]*5
        )

    def testMultiSavers(self):
        self.simulateAndCompare(*self._testMultiSaversParams())

    def testMultiSaversDB(self):
        dbPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}testMultiSavers.db')
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testMultiSaversParams(), dbPath, dbdct)

    def testFileBuffering(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        path = os.path.join(tempfile.gettempdir(), 'NewData.dat')

        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)
        rs = ResultSelector(self.newSim)
        # memory saving
        memsaver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count

        # file saving
        filesaver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
        filesaver.toFile(path, buffering=10)

        self.newSim.toSave(memsaver, filesaver)

        for rid in range(self.nbRuns):
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            for t in timepoints:
                self.newSim.run(t)
                memsaver.save()
                filesaver.save()
                if MPI._shouldWrite:
                    self.assertTrue((memsaver.data[-1]      == filesaver.data[-1]).all())
                    self.assertTrue((memsaver.data[rid]     == filesaver.data[rid]).all())
                    self.assertTrue((memsaver.data[-1, 0]   == filesaver.data[-1, 0]).all())
                    self.assertTrue((memsaver.data[rid, 0]  == filesaver.data[rid, 0]).all())
                    self.assertTrue((memsaver.data[-1, -1]  == filesaver.data[-1, -1]).all())
                    self.assertTrue((memsaver.data[rid, -1] == filesaver.data[rid, -1]).all())
                    self.assertTrue((memsaver.time[-1]      == filesaver.time[-1]).all())
                    self.assertTrue((memsaver.time[rid]     == filesaver.time[rid]).all())
                    self.assertEqual(memsaver.time[-1, 0]  , filesaver.time[-1, 0])
                    self.assertEqual(memsaver.time[rid, 0] , filesaver.time[rid, 0])
                    self.assertEqual(memsaver.time[-1, -1] , filesaver.time[-1, -1])
                    self.assertEqual(memsaver.time[rid, -1], filesaver.time[rid, -1])

                    self.assertTrue((memsaver.data[-1, :, 0] == filesaver.data[-1, :, 0]).all())
                    with np.testing.suppress_warnings() as sup:
                        sup.filter(np.VisibleDeprecationWarning)
                        self.assertTrue((memsaver.data[:, :, 0] == filesaver.data[:, :, 0]).all())
                    self.assertTrue((memsaver.data[-1,...] == filesaver.data[-1]).all())
                    self.assertTrue((memsaver.data[-1, ..., 0] == filesaver.data[-1, ..., 0]).all())
                    self.assertTrue((memsaver.data[-1, :, 0] == filesaver.data[-1, :, 0]).all())
                    self.assertTrue((memsaver.data[-1, 0, :] == filesaver.data[-1, 0, :]).all())
                    self.assertTrue((memsaver.data[-1, 0] == filesaver.data[-1, 0]).all())

                    with np.testing.suppress_warnings() as sup:
                        sup.filter(np.VisibleDeprecationWarning)
                        self.assertTrue((np.array(memsaver.data) == np.array(filesaver.data)).all())

                    with self.assertRaises(Exception):
                        memsaver.data[-1, :, 0, 0]
                    with self.assertRaises(Exception):
                        filesaver.data[-1, :, 0, 0]
                    with self.assertRaises(Exception):
                        filesaver.data[-1, :, 15]

                    with self.assertRaises(Exception):
                        filesaver.data[-1, ..., ...]

                    with self.assertRaises(Exception):
                        memsaver.time[-1, :, 0]
                    with self.assertRaises(Exception):
                        filesaver.time[-1, :, 0]

        if MPI._shouldWrite:
            self.assertSameData(memsaver.data[...], filesaver.data[...], [self.countRefVal]*3)

    def testResultSelectorLabels(self):

        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)
        rs = ResultSelector(self.newSim)

        ordSel = \
            rs.comp1.S1.Count <<\
            rs.comp1.S2.Conc <<\
            rs.LIST('comp1', 'comp2').S1.Count <<\
            rs.LIST(self.newGeom.comp1, self.newGeom.comp2).S2.Count <<\
            rs.comp1.CC.Count <<\
            rs.comp1.CC[self.newMdl.sus1, :].Count <<\
            rs.comp1.CC[..., self.newMdl.sus1].Count <<\
            rs.SUM(rs.comp1.LIST('S1', 'S2').Count) <<\
            rs.comp1.S1.Count + rs.comp2.S1.Count <<\
            (rs.comp1.LIST('S1', 'S2').Count + rs.comp2.LIST('S1', 'S2').Count) <<\
            rs.SUM(rs.comp1.LIST('S1', 'S2').Count * rs.comp2.LIST('S1', 'S2').Count)# <<\

        expOrdLbls = [
            'comp1.S1.Count',

            'comp1.S2.Conc',

            'comp1.S1.Count',
            'comp2.S1.Count',

            'comp1.S2.Count',
            'comp2.S2.Count',

            'comp1.CC.Count',
            'comp1.CC[sus1, :].Count',
            'comp1.CC[:, sus1].Count',

            'SUM(comp1.LIST(S1, S2).Count)',
            '(comp1.S1.Count + comp2.S1.Count)',

            '(comp1.S1.Count + comp2.S1.Count)',
            '(comp1.S2.Count + comp2.S2.Count)',

            'SUM((comp1.LIST(S1, S2).Count * comp2.LIST(S1, S2).Count))',
        ]

        lbls = ordSel.labels
        self.assertEqual(len(lbls), len(expOrdLbls))
        for l, el in zip(lbls, expOrdLbls):
            self.assertEqual(l, el)

        # Setting custom labels

        with self.assertRaises(Exception):
            ordSel.labels = ['test']

        with self.assertRaises(Exception):
            ordSel.labels = ['test'] * (len(expOrdLbls) + 1)

        newLabels = [f'label {i}' for i in range(len(expOrdLbls))]
        ordSel.labels = newLabels

        lbls = ordSel.labels
        self.assertEqual(len(lbls), len(newLabels))
        for l, el in zip(lbls, newLabels):
            self.assertEqual(l, el)


        # Operators
        ordSel = \
            (rs.comp1.LIST('S1', 'S2').Count + rs.comp2.LIST('S1', 'S2').Count) <<\
            (rs.comp1.LIST('S1', 'S2').Count - rs.comp2.LIST('S1', 'S2').Count) <<\
            (rs.comp1.LIST('S1', 'S2').Count * rs.comp2.LIST('S1', 'S2').Count) <<\
            (rs.comp1.LIST('S1', 'S2').Count / rs.comp2.LIST('S1', 'S2').Count) <<\
            (rs.comp1.LIST('S1', 'S2').Count ** rs.comp2.LIST('S1', 'S2').Count) <<\
            (123 + rs.comp2.LIST('S1', 'S2').Count) <<\
            (123 - rs.comp2.LIST('S1', 'S2').Count) <<\
            (123 * rs.comp2.LIST('S1', 'S2').Count) <<\
            (123 / rs.comp2.LIST('S1', 'S2').Count) <<\
            (rs.comp1.LIST('S1', 'S2').Count + 123) <<\
            (rs.comp1.LIST('S1', 'S2').Count - 123) <<\
            (rs.comp1.LIST('S1', 'S2').Count * 123) <<\
            (rs.comp1.LIST('S1', 'S2').Count / 123) <<\
            (rs.comp1.LIST('S1', 'S2').Count ** 123)

        expOrdLbls = []
        operators = ['+', '-', '*', '/', '**']
        tmp = ['(comp1.S1.Count {} comp2.S1.Count)', '(comp1.S2.Count {} comp2.S2.Count)']
        for op in operators:
            expOrdLbls += [t.format(op) for t in tmp]
        tmp = ['(123 {} comp2.S1.Count)', '(123 {} comp2.S2.Count)']
        for op in operators[:4]:
            expOrdLbls += [t.format(op) for t in tmp]
        tmp = ['(comp1.S1.Count {} 123)', '(comp1.S2.Count {} 123)']
        for op in operators:
            expOrdLbls += [t.format(op) for t in tmp]

        lbls = ordSel.labels
        self.assertEqual(len(lbls), len(expOrdLbls))
        for l, el in zip(lbls, expOrdLbls):
            self.assertEqual(l, el)

        # Unordered paths

        unordSel = \
            rs.ALL(Compartment, Patch).ALL(Species, Complex).Count

        expUnordLbls = [
            'comp1.S1.Count',
            'comp1.S2.Count',
            'comp1.CC.Count',
            'comp2.S1.Count',
            'comp2.S2.Count',
            'patch.Ex.Count',
            'patch.ExS1.Count',
            'patch.ExS2.Count',
            'patch.ExS1S2.Count',
        ]

        self.assertEqual(set(unordSel.labels), set(expUnordLbls))

        ##########
        unordSel = \
            rs.MATCH('^(co[mh]p.+)?(p[a-f]tch)?$').MATCH('^(Ex+)?(S[0-9])*$').Count

        expUnordLbls = [
            'comp1.S1.Count',
            'comp1.S2.Count',
            'comp2.S1.Count',
            'comp2.S2.Count',
            'patch.Ex.Count',
            'patch.ExS1.Count',
            'patch.ExS2.Count',
            'patch.ExS1S2.Count',
        ]

        self.assertEqual(set(unordSel.labels), set(expUnordLbls))

        ##########
        unordSel = \
            rs.comp1.LIST(*self.newMdl.CC).Count

        expUnordLbls = [f'comp1.{state}.Count' for state in self.newMdl.CC[...]]

        self.assertEqual(set(unordSel.labels), set(expUnordLbls))

        ##### Setting custom labels after start of simulation

        self.newSim.toSave(unordSel, dt = 0.1)

        self.newSim.newRun()

        with self.assertRaises(Exception):
            unordSel.labels = expUnordLbls

    
    def testMetaData(self, usingMesh=False):
        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)
        rs = ResultSelector(self.newSim)

        saver1 = rs.LIST('comp1', 'comp2').LIST('S1', 'S2').Count

        with self.assertRaises(TypeError):
            saver1.metaData[5] = [1, 2, 3, 4]
        with self.assertRaises(TypeError):
            saver1.metaData['test'] = 1
        with self.assertRaises(Exception):
            saver1.metaData['test'] = [1, 2]

        saver1.metaData['test'] = [1, 2, 3, 4]

        with self.assertRaises(KeyError):
            saver1.metaData['test2']

        self.assertEqual(list(saver1.metaData['test']), [1, 2, 3, 4])

        saver2 = rs.LIST('comp2', 'comp1').S1.Count
        saver2.metaData['test'] = [5, 6]

        saver12 = saver1 << saver2
        self.assertEqual(list(saver12.metaData['test']), [1, 2, 3, 4, 5, 6])

        saver2.metaData['test2'] = [7, 8]
        saver12 = saver1 << saver2
        self.assertEqual(list(saver12.metaData['test']), [1, 2, 3, 4, 5, 6])
        self.assertEqual(list(saver12.metaData['test2']), [None, None, None, None, 7, 8])

        # Automatic metadata
        self.assertEqual(list(saver1.metaData['loc_type']), [Compartment._locStr] * 4)
        self.assertEqual(list(saver1.metaData['loc_id']), ['comp1', 'comp1', 'comp2', 'comp2'])

        with self.assertWarns(Warning):
            saver1.metaData['loc_type'] = ['cmp'] * 4

        if usingMesh:
            n = len(self.newGeom.tets)
            saver3 = rs.TETS(self.newGeom.tets).S1.Count
            self.assertEqual(list(saver3.metaData['loc_type']), [TetReference._locStr] * n)
            self.assertEqual(list(saver3.metaData['loc_id']), list(range(n)))
            self.assertEqual(list(saver3.metaData['parent_loc_type']), [Compartment._locStr] * n)
            self.assertEqual(list(saver3.metaData['parent_loc_id']), [tet._getPhysicalLocation().name for tet in self.newGeom.tets])


    def testDBSaving(self):
        dbPath = os.path.join(tempfile.gettempdir(), f'{self.__class__.__name__}testDBSaving.db')
        if os.path.exists(dbPath) and MPI._shouldWrite:
            os.remove(dbPath)

        # self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)
        rs = ResultSelector(self.newSim)
        saver1 = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
        saver1 <<= rs.comp1.CC.sus2.Count
        saver1 <<= rs.SUM(rs.LIST('comp1', 'comp2').LIST('S1', 'S2').Count)

        saver1Mtdt = {'test': ['s1', 's2', 'exs1s2', 'CCsus2', 'sum'], 'test2':[1, 2, 3, 4, 5]}
        for key, lst in saver1Mtdt.items():
            saver1.metaData[key] = lst

        saver2 = rs.LIST('comp1', 'comp2').LIST('S1', 'S2').Count
        saver2lbls = ['c1S1', 'c1s2', 'c2s1', 'c2s2']
        saver2.labels = saver2lbls

        self.newSim.toSave(saver1, saver2, dt=self.deltaT)

        with SQLiteDBHandler(dbPath) as dbh:
            self.newSim.toDB(dbh, dbPath, val1=1, val2=2)
            self.newSim.newRun()
            self.newSim.run(self.shortEndTime)

        if MPI._shouldWrite:
            # Only rank 0 reads from the db, other ranks should not raise exceptiosn
            with self.assertRaises(Exception):
                self.newSim.toDB(dbh, dbPath, val1=1, val2=2)

        with SQLiteDBHandler(dbPath) as dbh:
            if MPI._shouldWrite:
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, dbPath, val1=3, val2=4)
            self.newSim.toDB(dbh, dbPath, val1=1, val2=2)
            self.newSim.newRun()
            self.newSim.run(self.shortEndTime)

        with SQLiteDBHandler(dbPath) as dbh:
            if MPI._shouldWrite:
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, f'{dbPath}_2', val3=[1, 2, 3])
            self.newSim.toDB(dbh, f'{dbPath}_2', val1=3, val2=4)
            self.newSim.newRun()
            self.newSim.run(self.shortEndTime)

        with SQLiteDBHandler(dbPath) as dbh:
            if MPI._shouldWrite:
                self.assertEqual(len(list(dbh)), 2)
                with self.assertRaises(TypeError):
                    dbh[42]
                with self.assertRaises(KeyError):
                    dbh['test']
                self.assertEqual(dbh[dbPath].parameters, {'val1':1, 'val2':2})
                self.assertEqual(dbh[f'{dbPath}_2'].parameters, {'val1':3, 'val2':4})
                self.assertEqual(len(dbh[dbPath].results), 2)
                self.assertEqual(dbh[dbPath].val1, 1)
                self.assertEqual(dbh[dbPath].val2, 2)
                with self.assertRaises(AttributeError):
                    dbh[dbPath].test
                self.assertEqual(dbh[dbPath].name, dbPath)

                sv1, sv2 = dbh[dbPath].results
                self.assertEqual(len(sv1.time), 2)
                self.assertEqual(len(sv2.time), 2)
                # Check labels
                self.assertEqual(saver1.labels, sv1.labels)
                self.assertEqual(saver2lbls, sv2.labels)
                # Check metadata
                saver1AutoMtdt = {
                    'loc_type': [Compartment._locStr] * 2 + [Patch._locStr] + [Compartment._locStr] * 2,
                    'loc_id': ['comp1', 'comp1', 'patch', 'comp1', None]
                }
                self.assertEqual({**saver1Mtdt, **saver1AutoMtdt}, sv1.metaData)
                self.assertEqual(saver2.metaData._dict, sv2.metaData)
            else:
                # Check that db access outside of rank 0 raises exceptions
                with self.assertRaises(Exception):
                    list(dbh)
                with self.assertRaises(Exception):
                    dbh[dbPath]

        # Manipulating resultselectors
        saver3 = rs.patch.LIST('ExS1', 'ExS2').Count
        self.newSim._resultSelectors = []
        self.newSim._nextSave = None
        self.newSim.toSave(saver1, saver3, dt=self.deltaT)
        with SQLiteDBHandler(dbPath) as dbh:
            if MPI._shouldWrite:
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, dbPath, val1=1, val2=2)

        self.newSim._resultSelectors = []
        self.newSim._nextSave = None
        self.newSim.toSave(saver1, saver2, saver3, dt=self.deltaT)
        with SQLiteDBHandler(dbPath) as dbh:
            if MPI._shouldWrite:
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, dbPath, val1=1, val2=2)

class TetSimDataSaving(test_model.TetTestModelFramework, SimDataSaving):
    """Test data access, setting, and saving with tetmeshes."""

    def setUp(self):
        test_model.TetTestModelFramework.setUp(self)
        SimDataSaving.setUp(self, False)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Tetexact', nmdl, ngeom, nrng, True)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        orng=srng.create('mt19937',512)
        orng.initialize(self.seed)
        osim = ssolver.Tetexact(omdl, ogeom, orng, True)
        osim.setEfieldDT(self.efielddt)
        return osim

    def testSimPathSetSyntax(self):
        super().testSimPathSetSyntax(usingMesh=True)

    def testSimPathGetSyntax(self, init=True, usingMesh=True):
        super().testSimPathGetSyntax(init=init, usingMesh=True)

    def _testDeltaTSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            c1 = 0
            for t in self.c1tets:
                c1 += sim.getTetCount(t.idx, 'S1')
            res.append(c1)
            p1 = 0
            for t in self.ptris:
                p1 += sim.getTriCount(t.idx, 'Ex')
            res.append(p1 * 2)
            return res

        def newSave(rs):
            saver = rs.SUM(rs.TETS(self.newGeom.comp1.tets).S1.Count) / 1.0
            tris = self.newGeom.patch.tris
            saver <<= rs.SUM(rs.TRIS(tris[0:len(tris)//2]).Ex.Count << rs.TRIS(tris[len(tris)//2:]).Ex.Count) * 2
            self.newSim.toSave(saver, dt=self.deltaT)
            return [saver]

        def newGetData(savers):
            saver = savers[0]
            saver.time[0]
            np.array(saver.time)
            saver.data[0,:]
            np.array(saver.data)
            return saver.data[:]

        def newRun(sim, rs):
            sim.run(self.endTime)

        explbls = [[None, None]]

        return (
            (timepoints, oldSave),
            (newSave, newGetData, newRun, explbls),
            [self.countRefVal]*2
        )

    def testMetaData(self):
        super().testMetaData(usingMesh=True)

    def testResultSelectorLabels(self):
        pass


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(SimDataSaving, "test"))
    all_tests.append(unittest.makeSuite(TetSimDataSaving, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

