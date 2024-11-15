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

""" Unit tests for data saving."""

import functools
import importlib.util
import itertools
import numpy as np
import operator
import os
import re
import shutil
import sys
import tempfile
import threading
import time
import unittest

import steps
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

from . import base_model


class ResultSelectorTests(base_model.TestModelFramework):
    """Test result selector functionality that does not depend on a specific solver"""

    def setUp(self, callParent=True):
        if callParent:
            super().setUp()

        self.newMdl = self.get_API2_Mdl()
        self.newGeom = self.get_API2_Geom(self.newMdl)
        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Wmdirect', nmdl, ngeom, nrng)
        return nsim

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

    def testResultSelectorLabels(self):
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

    def testDescription(self):
        rs = ResultSelector(self.newSim)

        saver1 = rs.LIST('comp1', 'comp2').LIST('S1', 'S2').Count

        self.assertEqual(saver1.description, 'LIST(comp1, comp2).LIST(S1, S2).Count')

        with self.assertRaises(TypeError):
            saver1.description = 34

        saver1.description = 'my custom description'
        self.assertEqual(saver1.description, 'my custom description')

        saverDescr = [
            (rs.comp1.S1.Conc, 'comp1.S1.Conc'),
            (rs.ALL(Compartment).S1.Count, 'ALL(Compartment).S1.Count'),
            (rs.MATCH(r'comp\d').S1.Count, r'MATCH(comp\d).S1.Count'),
            (rs.LIST(self.newGeom.comp1, self.newGeom.comp2).S2.Count, 'LIST(comp1, comp2).S2.Count'),
            (rs.comp1.S1.Count << rs.comp1.S2.Count, 'comp1.S1.Count, comp1.S2.Count'),
            (rs.comp1.S1.Count + rs.comp1.S2.Count, '(comp1.S1.Count + comp1.S2.Count)'),
            (rs.comp1.S1.Count - rs.comp1.S2.Count, '(comp1.S1.Count - comp1.S2.Count)'),
            (rs.comp1.S1.Count * rs.comp1.S2.Count, '(comp1.S1.Count * comp1.S2.Count)'),
            (rs.comp1.S1.Count / rs.comp1.S2.Count, '(comp1.S1.Count / comp1.S2.Count)'),
            (rs.comp1.S1.Count ** rs.comp1.S2.Count, '(comp1.S1.Count ** comp1.S2.Count)'),
            (10 * rs.comp1.S1.Count + rs.comp1.S2.Count, '((10 * comp1.S1.Count) + comp1.S2.Count)'),
            (rs.SUM(rs.comp1.ALL(Species).Count), 'SUM(comp1.ALL(Species).Count)'),
            (rs.MAX(rs.comp1.ALL(Species).Count), 'MAX(comp1.ALL(Species).Count)'),
            (rs.MIN(rs.comp1.ALL(Species).Count), 'MIN(comp1.ALL(Species).Count)'),
        ]

        for saver, descr in saverDescr:
            self.assertEqual(saver.description, descr)

    def _checkFeature(self, objStr, obj, tpe, version, firstAdded, noneok=False, exceptTpe=NotImplementedError):
        if version >= firstAdded:
            values = eval(objStr)
            self.assertGreater(len(values), 0)
            self.assertTrue(all(isinstance(v, tpe) or (v is None and noneok) for v in values))
        else:
            with self.assertRaises(exceptTpe):
                eval(objStr)

    def _testLoadingResultSelector(self, rs, fmt, version):
        floatTpe = (float, np.float32, np.float64)
        self._checkFeature("obj.time[0]", rs, floatTpe, version, (3,6,0))
        self._checkFeature("obj.data[0,:,0]", rs, floatTpe, version, (3,6,0))
        self._checkFeature("obj.labels", rs, str, version, (3,6,0))
        self._checkFeature("obj.metaData['loc_type']", rs, str, version, (4,0,0), noneok=True)
        self._checkFeature("obj.description", rs, str, version, (5,0,1))

    def _testLoadingDBGroup(self, db, fmt, version):
        group = db.get()

        if fmt == 'h5':
            self.assertEqual(version, group._version)

        self._checkFeature("obj.name", group, str, version, (4,0,0))
        self._checkFeature("obj.results", group, steps.saving._ReadOnlyResultSelector, version, (4,0,0))
        self._checkFeature("obj.parameters", group, str, version, (4,0,0))
        # Staticdata not implemented for SQLite db
        self._checkFeature("obj.staticData.keys()", group, str, version if fmt == 'h5' else (0, 0, 0), (5,0,0),
                           exceptTpe=(NotImplementedError, steps.saving.UnavailableDataError))

        rs, *_ = group.results
        self._testLoadingResultSelector(rs, fmt, version)

    def testLoadOlderFiles(self):
        old_dirs = os.path.join(base_model.FILEDIR, '..', 'saved_files')
        for dir in os.listdir(old_dirs):
            dirpath = os.path.join(old_dirs, dir)
            if os.path.isdir(dirpath):
                try:
                    version = steps.utils.Versioned._parseVersion(dir)
                except ValueError:
                    continue
                for file in os.listdir(dirpath):
                    *_, ext = file.split('.')
                    path = os.path.join(dirpath, file)
                    with self.subTest(fileFormat=ext, fileVersion=version):
                        if ext == 'dat':
                            rs = ResultSelector.FromFile(path)
                            self._testLoadingResultSelector(rs, ext, version)
                        elif ext == 'db':
                            with SQLiteDBHandler(path) as db:
                                self._testLoadingDBGroup(db, ext, version)
                        elif ext == 'h5':
                            with HDF5Handler(path[:-3]) as db:
                                self._testLoadingDBGroup(db, ext, version)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testHDF5MultiReaderBrokenFiles(self):
        saved_dir = os.path.join(base_model.FILEDIR, '..', 'saved_files')
        broken_dir = os.path.join(saved_dir, 'broken')

        # All files are broken, nothing to load so we should raise FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            with self.assertWarns(RuntimeWarning):
                hdf = HDF5MultiFileReader(broken_dir, recursive=True)

        # Saved files folder contains some loadable files, we should only raise warnings
        with self.assertWarns(RuntimeWarning):
            hdf = HDF5MultiFileReader(saved_dir, recursive=True)


class SimDataSaving(base_model.TestModelFramework):
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

        self.createdFiles = set()
        self.createdDirectories = set()

        self.sqliteArgs = []
        self.sqliteKwargs = dict(commitFreq=100)

        self.hdf5Args = []
        self.hdf5Kwargs = {}

        self.xdmfArgs = []
        self.xdmfKwargs = {}

    def tearDown(self):
        super().tearDown()
        for path in self.createdFiles:
            if os.path.isfile(path):
                os.remove(path)
        for path in self.createdDirectories:
            # Delete directories, starting with the most nested
            if os.path.isdir(path):
                shutil.rmtree(path)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Wmdirect', nmdl, ngeom, nrng)
        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        orng=srng.create('mt19937',512)
        orng.initialize(self.seed)
        osim = ssolver.Wmdirect(omdl, ogeom, orng)
        return osim

    def testSimPathGetSyntax(self, init=True):
        self._test_API2_SimPathGetSyntax(self.newSim, init=init)

    def testSimPathSetSyntax(self):
        self._test_API2_SimPathSetSyntax(self.newSim)

    def simulateAndCompare(self, oldSimToSave, newSimToSave, refVal, dbCls=None, dbPath=None, dbUID=None,
            dbArgs=None, dbKwargs=None, dbGroupArgs={}):
        if self.oldSim is not None:
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

        if dbCls is None:
            for rid in range(self.nbRuns):
                self.newSim.newRun()
                self.init_API2_sim(self.newSim)
                runFunc(self.newSim, savers)
            if MPI._shouldWrite:
                newDat = getDataFunc(savers)
                labelList = [s.labels for s in savers]
                descriptions = [s.description for s in savers]
        else:
            # Test database saving
            with dbCls(dbPath, *dbArgs, **dbKwargs) as dbh:
                group = self.newSim.toDB(dbh, dbUID, **dbGroupArgs)
                if isinstance(dbh, HDF5Handler) and MPI._shouldWrite:
                    group.staticData['staticDataTest'] = [1, 2, 3]
                for rid in range(self.nbRuns):
                    self.newSim.newRun()
                    self.init_API2_sim(self.newSim)
                    runFunc(self.newSim, savers)

                if MPI._shouldWrite:
                    newDat = getDataFunc(savers)
                    labelList = [s.labels for s in savers]
                    descriptions = [s.description for s in savers]

                self.createdFiles |= set(dbh._getFilePaths())

            # Data checking is only done in rank 0
            if MPI._shouldWrite:
                with dbCls(dbPath, *dbArgs, **dbKwargs) as dbh:
                    for name, val in dbGroupArgs.items():
                        self.assertEqual(getattr(dbh[dbUID], name), val)
                    # parameter selection
                    self.assertEqual(dbh.parameters.keys(), dbGroupArgs.keys())
                    for name, val in dbh.parameters.items():
                        self.assertEqual(len(val), 1)
                        self.assertEqual(list(val)[0], dbGroupArgs[name])
                    groups = dbh.filter(**dbGroupArgs)
                    self.assertEqual(len(groups), 1)
                    self.assertEqual(list(groups)[0].name, dbUID)
                    # results
                    loadedSavers = dbh[dbUID].results
                    newDat2 = getDataFunc(loadedSavers)
                    labelList2 = [s.labels for s in loadedSavers]
                    descriptions2 = [s.description for s in loadedSavers]

            # Check that the loaded data is identical to the one extracted during simulation
            if MPI._shouldWrite:
                self.assertSameData(newDat, newDat2, refVal)
                self.assertEqual(descriptions, descriptions2)
                for llst, explbl in zip(labelList, labelList2):
                    self.assertEqual(len(llst), len(explbl))
                    for rslbl, elbl in zip(llst, explbl):
                        self.assertEqual(rslbl, elbl)

        if MPI._shouldWrite:
            if self.oldSim is not None:
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
            res.append(sim.getCompSpecCount('comp1', 'S1'))
            res.append(sim.getCompSpecCount('comp1', 'S2'))
            res.append(sim.getPatchSpecCount('patch', 'ExS1S2'))
            res.append(sim.getCompSpecCount('comp1', 'CCsus12'))
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

        explbls = [['comp1.S1.Count', 'comp1.S2.Count', 'patch.ExS1S2.Count', 'comp1.CC[sus1, sus2].Count']]

        return (
            (timepoints, oldSave),
            (newSave, newGetData, newRun, explbls),
            [self.countRefVal]*3,
        )

    def testDeltaTSaving(self):
        self.simulateAndCompare(*self._testDeltaTSavingParams())

    def _testDeltaTSavingDB(self, dbCls, dbArgs, dbKwargs, dbUID='GroupId', suffix=''):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testDeltaTSaving', suffix=suffix)
        self.createdFiles.add(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testDeltaTSavingParams(), dbCls, dbPath, dbUID, dbArgs, dbKwargs, dbdct)

    def testDelataTSavingSQLite(self):
        self._testDeltaTSavingDB(SQLiteDBHandler, self.sqliteArgs, self.sqliteKwargs, suffix='.db')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testDeltaTSavingHDF5(self):
        self._testDeltaTSavingDB(HDF5Handler, self.hdf5Args, self.hdf5Kwargs)

    def _testTimepointSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompSpecCount('comp1', 'S1'))
            res.append(sim.getCompSpecCount('comp1', 'S2'))
            res.append(sim.getPatchSpecCount('patch', 'ExS1S2'))
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

    def _testTimepointSavingDB(self, dbCls, dbArgs, dbKwargs, dbUID='GroupId', suffix=''):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testTimepointSaving', suffix='.db')
        self.createdFiles.add(dbPath)
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testTimepointSavingParams(), dbCls, dbPath, dbUID, dbArgs, dbKwargs, dbdct)

    def testTimepointSavingSQLite(self):
        self._testTimepointSavingDB(SQLiteDBHandler, self.sqliteArgs, self.sqliteKwargs, suffix='.db')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testTimepointSavingHDF5(self):
        self._testTimepointSavingDB(HDF5Handler, self.hdf5Args, self.hdf5Kwargs)

    def _testUnspecifiedSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompSpecCount('comp1', 'S1'))
            res.append(sim.getCompSpecCount('comp1', 'S2'))
            res.append(sim.getPatchSpecCount('patch', 'ExS1S2'))
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

    def _testUnspecifiedSavingDB(self, dbCls, dbArgs, dbKwargs, dbUID='GroupId', suffix=''):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testUnspecifiedSaving', suffix='.db')
        self.createdFiles.add(dbPath)
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testUnspecifiedSavingParams(), dbCls, dbPath, dbUID, dbArgs, dbKwargs, dbdct)

    def testUnspecifiedSavingSQLite(self):
        self._testUnspecifiedSavingDB(SQLiteDBHandler, self.sqliteArgs, self.sqliteKwargs, suffix='.db')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testUnspecifiedSavingHDF5(self):
        self._testUnspecifiedSavingDB(HDF5Handler, self.hdf5Args, self.hdf5Kwargs)

    def _testFileSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        _, path = tempfile.mkstemp(prefix=f'NewData', suffix='.dat')
        self.createdFiles.add(path)
        def oldSave(sim):
            res = []
            res.append(sim.getCompSpecCount('comp1', 'S1'))
            res.append(sim.getCompSpecCount('comp1', 'S2'))
            res.append(sim.getPatchSpecCount('patch', 'ExS1S2'))
            res.append(sim.getCompSpecCount('comp1', 'CCsus12') + 2*sim.getCompSpecCount('comp1', 'CCsus12'))
            tmp = sim.getCompSpecCount('comp1', 'S1') + sim.getCompSpecCount('comp1', 'S2') + sim.getCompSpecCount('comp2', 'S1') + sim.getCompSpecCount('comp2', 'S2')
            tmp /= (sim.getPatchSpecCount('patch', 'ExS1S2') + sim.getPatchSpecCount('patch', 'Ex'))
            res.append(tmp)
            return res

        def newSave(rs):
            saver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count
            saver <<= rs.comp1.CC.sus2.Count
            saver <<= rs.SUM(rs.LIST('comp1', 'comp2').LIST('S1', 'S2').Count) / rs.SUM(rs.patch.ExS1S2.Count << rs.patch.Ex.Count)
            self.newSim.toSave(saver, dt=self.deltaT)
            saver.toFile(path)
            return [saver]

        def newGetData(savers, fromFile=True):
            saver = savers[0]
            saver.time[0]
            np.array(saver.time)
            saver.data[0,:]
            np.array(saver.data)
            if fromFile:
                rs = ResultSelector.FromFile(path)
            else:
                rs = saver
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

    def testFileSavingOlderVersions(self):
        version = steps.__version__
        steps.__version__ = '3.5.0'
        self.simulateAndCompare(*self._testFileSavingParams())
        steps.__version__ = version

    def _testFileSavingDB(self, dbCls, dbArgs, dbKwargs, dbUID='GroupId', suffix=''):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testFileSaving', suffix='.db')
        self.createdFiles.add(dbPath)
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        params = list(self._testFileSavingParams())
        ngd = params[1][1]
        params[1] = params[1][0:1] + (lambda savers: ngd(savers, fromFile=False),) + params[1][2:]
        self.simulateAndCompare(*params, dbCls, dbPath, dbUID, dbArgs, dbKwargs, dbdct)

    def testFileSavingSQLite(self):
        self._testFileSavingDB(SQLiteDBHandler, self.sqliteArgs, self.sqliteKwargs, suffix='.db')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testFileSavingHDF5(self):
        self._testFileSavingDB(HDF5Handler, self.hdf5Args, self.hdf5Kwargs)

    def testPartialDataReading(self):
        _, path = tempfile.mkstemp(prefix=f'NewData', suffix='.dat')
        self.createdFiles.add(path)

        # Data reading function
        exitThread = False
        def readData(path):
            while not exitThread:
                rs = ResultSelector.FromFile(path)
                nbIndexErr = 0
                try:
                    rs.data[self.nbRuns-1, 0]
                except Exception:
                    nbIndexErr += 1
                try:
                    rs.time[self.nbRuns-1, 0]
                except Exception:
                    nbIndexErr += 1

                try:
                    n = len(rs.data) - 1
                except Exception:
                    continue
                if n < self.nbRuns - 1:
                    self.assertEqual(nbIndexErr, 2)
                    try:
                        rs.data[n, -1]
                        rs.time[n, -1]
                    except Exception:
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
            is_alive = thr.is_alive()
            if is_alive:
                exitThread = True
                thr.join(1)
            self.assertFalse(is_alive)

    def _testMultiSaversParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            res.append(sim.getCompSpecCount('comp1', 'S1'))
            res.append(sim.getCompSpecCount('comp1', 'S2'))
            res.append(sim.getPatchSpecCount('patch', 'ExS1S2'))
            res.append(sim.getCompSpecCount('comp2', 'S1'))
            res.append(sim.getCompSpecCount('comp2', 'S2'))
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

    def _testMultiSaversDB(self, dbCls, dbArgs, dbKwargs, dbUID='GroupId', suffix=''):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testMultiSavers', suffix='.db')
        self.createdFiles.add(dbPath)
        if MPI._shouldWrite and os.path.exists(dbPath):
            os.remove(dbPath)
        dbdct = {'test1': 'txtvalue', 'test2': 42, 'test3': 5.245}
        self.simulateAndCompare(*self._testMultiSaversParams(), dbCls, dbPath, dbUID, dbArgs, dbKwargs, dbdct)

    def testMultiSaversSQLite(self):
        self._testMultiSaversDB(SQLiteDBHandler, self.sqliteArgs, self.sqliteKwargs, suffix='.db')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testMultiSaversHDF5(self):
        self._testMultiSaversDB(HDF5Handler, self.hdf5Args, self.hdf5Kwargs)

    def testFileBuffering(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        _, path = tempfile.mkstemp(prefix=f'NewData', suffix='.dat')
        self.createdFiles.add(path)

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
                        sup.filter(np.exceptions.VisibleDeprecationWarning)
                        self.assertTrue((memsaver.data[:, :, 0] == filesaver.data[:, :, 0]).all())
                    self.assertTrue((memsaver.data[-1,...] == filesaver.data[-1]).all())
                    self.assertTrue((memsaver.data[-1, ..., 0] == filesaver.data[-1, ..., 0]).all())
                    self.assertTrue((memsaver.data[-1, :, 0] == filesaver.data[-1, :, 0]).all())
                    self.assertTrue((memsaver.data[-1, 0, :] == filesaver.data[-1, 0, :]).all())
                    self.assertTrue((memsaver.data[-1, 0] == filesaver.data[-1, 0]).all())

                    with np.testing.suppress_warnings() as sup:
                        sup.filter(np.exceptions.VisibleDeprecationWarning)
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

    def testMetaData(self):
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

        if self.useMesh:
            n = len(self.newGeom.tets)
            saver3 = rs.TETS(self.newGeom.tets).S1.Count
            self.assertEqual(list(saver3.metaData['loc_type']), [TetReference._locStr] * n)
            self.assertEqual(list(saver3.metaData['loc_id']), list(range(n)))
            self.assertEqual(list(saver3.metaData['parent_loc_type']), [Compartment._locStr] * n)
            self.assertEqual(list(saver3.metaData['parent_loc_id']), [tet._getPhysicalLocation().name for tet in self.newGeom.tets])

    def testSQLiteDBSaving(self):
        self._testDBSaving(SQLiteDBHandler, 'SQLite', '.db')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testHDF5Saving(self):
        self._testDBSaving(HDF5Handler, 'HDF5')

    def _testDBSaving(self, dbCls, dbName, dbUID='GroupID', suffix=''):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}test{dbName}DBSaving', suffix=suffix)
        self.createdFiles.add(dbPath)

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

        with dbCls(dbPath) as dbh:
            # No slashes in group uid
            with self.assertRaises(ValueError):
                self.newSim.toDB(dbh, f'{dbUID}/test')

            # Default group uid
            self.newSim.toDB(dbh, val1=4, val2=5)
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

            # Specific group uid
            group = self.newSim.toDB(dbh, dbUID, val1=1, val2=2)

            # Static data
            if MPI._shouldWrite:
                # Also check that parameters can be accessed
                self.assertCountEqual(group.parameters.items(), dict(val1=1, val2=2).items())
                self.assertEqual(group.name, dbUID)
                if dbCls == SQLiteDBHandler:
                    with self.assertRaises(NotImplementedError):
                        group.staticData
                else:
                    group.staticData['teststatic1'] = list(range(10))
                    group.staticData['teststatic1'] = list(range(10))
                    with self.assertRaises(Exception):
                        group.staticData['teststatic1'] = list(range(5))
                    with self.assertRaises(KeyError):
                        group.staticData[42] = [1, 2, 3]
                    with self.assertRaises(TypeError):
                        group.staticData['teststatic2'] = set([1, 2, 3])
            else:
                self.assertIsNone(group)

            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

        if MPI._shouldWrite:
            # Only rank 0 reads from the db, other ranks should not raise exceptiosn
            with self.assertRaises(Exception):
                self.newSim.toDB(dbh, dbUID, val1=1, val2=2)

        with dbCls(dbPath) as dbh:
            if MPI._shouldWrite:
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, dbUID, val1=3, val2=4)
            group = self.newSim.toDB(dbh, dbUID, val1=1, val2=2)

            # Static data
            if MPI._shouldWrite:
                # Also check that parameters can be accessed
                self.assertCountEqual(group.parameters.items(), dict(val1=1, val2=2).items())
                self.assertEqual(group.name, dbUID)
                if dbCls == SQLiteDBHandler:
                    with self.assertRaises(NotImplementedError):
                        group.staticData
                else:
                    group.staticData['teststatic2'] = {i: 'a'*i for i in range(5)}
                    with self.assertRaises(Exception):
                        group.staticData['teststatic1'] = list(range(5))
                    with self.assertRaises(KeyError):
                        group.staticData[42] = [1, 2, 3]
                    with self.assertRaises(TypeError):
                        group.staticData['teststatic3'] = set([1, 2, 3])
            else:
                self.assertIsNone(group)

            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

        with dbCls(dbPath) as dbh:
            if MPI._shouldWrite:
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, f'{dbUID}_2', val3=[[1, 5], 2, 3])
            self.newSim.toDB(dbh, f'{dbUID}_2', val1=3, val2=4)
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

        # Read data
        with dbCls(dbPath) as dbh:
            if MPI._shouldWrite:
                self.assertEqual(len(dbh), 3)
                with self.assertRaises(TypeError):
                    dbh[42]
                with self.assertRaises(KeyError):
                    dbh['test']
                # dbh properties
                defaultName = DatabaseHandler._DEFAULT_GROUP_NAME.format(0)
                self.assertEqual(dbh.parameters, dict(val1=set([1, 3, 4]), val2=set([2, 4, 5])))
                self.assertEqual(dbh.filter(val1=1, val2=2), set([dbh[dbUID]]))
                self.assertEqual(dbh.filter(val1=3, val2=4), set([dbh[f'{dbUID}_2']]))
                self.assertEqual(dbh.filter(val1=4, val2=5), set([dbh[defaultName]]))
                self.assertEqual(dbh.filter(val1=1, val2=4), set())
                self.assertEqual(dbh.filter(val1=5), set())
                with self.assertRaises(KeyError):
                    dbh.filter(val4=5)
                self.assertEqual(dbh.filter(), set(dbh))
                self.assertEqual(dbh.get(val1=1, val2=2), dbh[dbUID])
                with self.assertRaises(Exception):
                    dbh.get(val1=1, val2=4)
                # run group properties
                self.assertEqual(dbh[dbUID].parameters, {'val1':1, 'val2':2})
                self.assertEqual(dbh[f'{dbUID}_2'].parameters, {'val1':3, 'val2':4})
                self.assertEqual(dbh[defaultName].parameters, {'val1':4, 'val2':5})
                self.assertEqual(len(dbh[dbUID].results), 2)
                self.assertEqual(dbh[dbUID].val1, 1)
                self.assertEqual(dbh[dbUID].val2, 2)
                with self.assertRaises(AttributeError):
                    dbh[dbUID].test
                self.assertEqual(dbh[dbUID].name, dbUID)
                self.assertEqual(dbh[defaultName].name, defaultName)
                # Static data
                if dbCls != SQLiteDBHandler:
                    self.assertEqual(dbh[dbUID].staticData['teststatic1'], list(range(10)))
                    self.assertEqual(dbh[dbUID].staticData['teststatic2'], {i: 'a'*i for i in range(5)})

                sv1, sv2 = dbh[dbUID].results
                self.assertEqual(len(sv1.time), 2)
                self.assertEqual(len(sv2.time), 2)
                # Check labels
                self.assertEqual(saver1.labels, sv1.labels)
                self.assertEqual(saver2lbls, sv2.labels)
                # Check metadata
                saver1AutoMtdt = {
                    'loc_type': [Compartment._locStr] * 2 + [Patch._locStr] + [Compartment._locStr] * 2,
                    'loc_id': ['comp1', 'comp1', 'patch', 'comp1', None],
                    'obj_type': ['Spec', 'Spec', 'Spec', '', 'Spec'],
                    'obj_id': ['S1', 'S2', 'ExS1S2', 'sus2', None],
                    'property': ['Count', 'Count', 'Count', 'Count', 'Count'],
                }
                self.assertEqual({**saver1Mtdt, **saver1AutoMtdt}, sv1.metaData)
                self.assertEqual(saver2.metaData._dict, sv2.metaData)
            else:
                # Check that db access outside of rank 0 raises exceptions
                with self.assertRaises(Exception):
                    list(dbh)
                with self.assertRaises(Exception):
                    dbh[dbUID]

        # Manipulating resultselectors
        saver3 = rs.patch.LIST('ExS1', 'ExS2').Count
        self.newSim._resultSelectors = []
        self.newSim._nextSave = None
        self.newSim.toSave(saver1, saver3, dt=self.deltaT)
        with dbCls(dbPath) as dbh:
            if MPI._shouldWrite or (isinstance(dbh, HDF5Handler) and self.newSim._isDistributed()):
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            else:
                self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            self.createdFiles |= set(dbh._getFilePaths())

        self.newSim._resultSelectors = []
        self.newSim._nextSave = None
        self.newSim.toSave(saver1, saver2, saver3, dt=self.deltaT)
        with dbCls(dbPath) as dbh:
            if MPI._shouldWrite or (isinstance(dbh, HDF5Handler) and self.newSim._isDistributed()):
                with self.assertRaises(Exception):
                    self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            else:
                self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            self.createdFiles |= set(dbh._getFilePaths())

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testXDMFWithoutMesh(self, dbUID='GroupID'):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testXDMFWithoutMeshSaving', suffix='')
        self.createdFiles.add(dbPath)

        # Check that XDMFHandler does not work with well-mixed simulations
        if not isinstance(self, TetSimDataSaving):
            with XDMFHandler(dbPath) as dbh:
                with self.assertRaises(TypeError):
                    self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
                self.createdFiles |= set(dbh._getFilePaths())

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testHDF5MultiReader(self):
        dirParams = {
            'param1': [1, 2, 3],
            'param2': ['a', 'b', 'c'],
            'param3': [0.1, 0.2, 0.3],
        }
        fileParams = {
            'param4': [4, 5, 6],
        }
        innerParams = {
            'param5': [0.1, 0.2]
        }
        nDir = functools.reduce(operator.mul, map(len, dirParams.values()), 1)
        nFile = functools.reduce(operator.mul, map(len, fileParams.values()), 1)
        nInner = functools.reduce(operator.mul, map(len, innerParams.values()), 1)

        dirPath = ''
        if MPI._shouldWrite:
            dirPath = tempfile.mkdtemp(prefix=f'{self.__class__.__name__}testHDF5MultiReader')
            self.createdDirectories.add(dirPath)
        if MPI._usingMPI and MPI.nhosts > 1:
            import mpi4py.MPI
            dirPath = mpi4py.MPI.COMM_WORLD.bcast(dirPath, root=0)
        if MPI._shouldWrite:
            emptyDir = tempfile.mkdtemp(prefix='empty', dir=dirPath)

        rs = ResultSelector(self.newSim)
        saver = rs.comp1.S1.Count << rs.comp1.S2.Count << rs.patch.ExS1S2.Count

        self.newSim.toSave(saver, dt=self.deltaT)

        allPaths = []
        for values in itertools.product(*dirParams.values()):
            kwargs = {name: v for name, v in zip(dirParams.keys(), values)}
            subDirPath = os.path.join(dirPath, '/'.join(f'{k}_{v}' for k, v in kwargs.items()))
            if MPI._shouldWrite:
                os.makedirs(subDirPath, exist_ok=True)
            if MPI._usingMPI and MPI.nhosts > 1:
                mpi4py.MPI.COMM_WORLD.barrier()
            # For each directory, add several files
            for values2 in itertools.product(*fileParams.values()):
                kwargs2 = {name: v for name, v in zip(fileParams.keys(), values2)}
                dbPath = os.path.join(subDirPath, '__'.join(f'{k}_{v}' for k, v in kwargs2.items()))
                allPaths.append(dbPath)

                with HDF5Handler(dbPath) as hdf:
                    # In each file, add several groups
                    for values3 in itertools.product(*innerParams.values()):
                        kwargs3 = {name: v for name, v in zip(innerParams.keys(), values3)}
                        self.newSim.toDB(hdf, **kwargs, **kwargs2, **kwargs3)
                        self.newSim.newRun()
                        self.init_API2_sim(self.newSim)
                        self.newSim.run(self.shortEndTime)
                        self.createdFiles |= set(hdf._getFilePaths())

        # Only read data in rank 0
        if MPI.rank == 0:
            # Loading from emptyDirectory
            with self.assertRaises(FileNotFoundError):
                hdf = HDF5MultiFileReader(emptyDir)

            # Non-existent directory
            with self.assertRaises(FileNotFoundError):
                hdf = HDF5MultiFileReader(emptyDir + '_2')

            # Empty list of files
            with self.assertRaises(FileNotFoundError):
                hdf = HDF5MultiFileReader([])

            allParams = {**dirParams, **fileParams, **innerParams}
            def checkAllFiles(hdf):
                self.assertEqual(len(hdf), nDir * nFile * nInner)
                self.assertEqual(hdf.parameters, {k: set(v) for k, v in allParams.items()})
                self.assertEqual(set(hdf.filter(param5=0.1)), set(hdf[DatabaseHandler._DEFAULT_GROUP_NAME.format(0)]))
                self.assertEqual(set(hdf.filter(param5=0.2)), set(hdf[DatabaseHandler._DEFAULT_GROUP_NAME.format(1)]))
                for group in hdf:
                    self.assertTrue(re.match(re.sub(r'{.*}', r'\\d+', DatabaseHandler._DEFAULT_GROUP_NAME), group.name) is not None)
                    self.assertTrue(all(param in allParams and value in allParams[param] for param, value in group.parameters.items()))
                self.assertEqual(len(hdf.filter(**{k: v[0] for k, v in dirParams.items()})), nFile * nInner)
                self.assertEqual(len(hdf.filter(**{k: v[0] for k, v in fileParams.items()})), nDir * nInner)
                self.assertEqual(len(hdf.filter(**{k: v[0] for k, v in innerParams.items()})), nDir * nFile)
                self.assertEqual(len(hdf.filter(param1=34)), 0)
                with self.assertRaises(KeyError):
                    hdf.filter(param7=34)
                with self.assertRaises(KeyError):
                    hdf['test42']
                with self.assertRaises(Exception):
                    hdf.get(param1=1)
                with self.assertRaises(Exception):
                    hdf.get(param1=17)
                with self.assertRaises(KeyError):
                    hdf.get(testparam=3)

            # Loading from list of file prefixes
            with HDF5MultiFileReader(allPaths) as hdf:
                checkAllFiles(hdf)

            # Loading recursively from top directory
            with HDF5MultiFileReader(dirPath, recursive=True) as hdf:
                checkAllFiles(hdf)

            # Loading from top directory (empty)
            with self.assertRaises(FileNotFoundError):
                hdf = HDF5MultiFileReader(dirPath)

            # Loading from sub directory
            with HDF5MultiFileReader(os.path.join(dirPath, 'param1_1'), recursive=True) as hdf:
                self.assertEqual(len(hdf), nDir * nFile * nInner // len(dirParams['param1']))

        if MPI._usingMPI and MPI.nhosts > 1:
            mpi4py.MPI.COMM_WORLD.barrier()

        # Add files with different parameters
        with HDF5Handler(os.path.join(dirPath, 'singlefile')) as hdf:
            self.newSim.toDB(hdf, 'singleUID', param6=3, param7=8)
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(hdf._getFilePaths())

        if MPI._shouldWrite:
            # Load only single file
            singleParams = {'param6': set([3]), 'param7': set([8])}
            with HDF5MultiFileReader(dirPath) as hdf:
                self.assertEqual(len(hdf), 1)
                self.assertEqual(hdf.parameters, singleParams)
                self.assertEqual(set(hdf.filter(param6=3, param7=8)), set(hdf))
                self.assertEqual(hdf['singleUID'], hdf.get())

            # Load everything
            with HDF5MultiFileReader(dirPath, recursive=True) as hdf:
                self.assertEqual(len(hdf), nDir * nFile * nInner + 1)
                self.assertEqual(hdf.parameters, {**{k: set(v) for k, v in allParams.items()}, **singleParams})
                self.assertEqual(hdf['singleUID'], hdf.get(param6=3))
                self.assertEqual(hdf['singleUID'], hdf.get(param6=3, param7=8))
                self.assertEqual(set([hdf['singleUID']]), hdf.filter(param6=3, param7=8))
                self.assertEqual(len(hdf.filter(**{k: v[0] for k, v in dirParams.items()})), nFile * nInner)

        if MPI._usingMPI and MPI.nhosts > 1:
            mpi4py.MPI.COMM_WORLD.barrier()


class TetSimDataSaving(base_model.TetTestModelFramework, SimDataSaving):
    """Test data access, setting, and saving with tetmeshes."""

    def setUp(self, callParent=True):
        if callParent:
            base_model.TetTestModelFramework.setUp(self)
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

    def _testDeltaTSavingParams(self):
        timepoints = np.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            c1 = 0
            for t in self.oc1tets:
                c1 += sim.getTetSpecCount(t, 'S1')
            res.append(c1)
            p1 = 0
            for t in self.optris:
                p1 += sim.getTriSpecCount(t, 'Ex')
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
        super().testMetaData()

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testDeltaTSavingXDMF(self, dbUID='GroupId1'):
        self._testDeltaTSavingDB(XDMFHandler, self.xdmfArgs, self.xdmfKwargs, dbUID=dbUID)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testTimepointSavingXDMF(self, dbUID='GroupId2'):
        self._testTimepointSavingDB(XDMFHandler, self.hdf5Args, self.hdf5Kwargs, dbUID=dbUID)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testUnspecifiedSavingXDMF(self, dbUID='GroupId3'):
        self._testUnspecifiedSavingDB(XDMFHandler, self.hdf5Args, self.hdf5Kwargs, dbUID=dbUID)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testFileSavingXDMF(self, dbUID='GroupId4'):
        self._testFileSavingDB(XDMFHandler, self.hdf5Args, self.hdf5Kwargs, dbUID=dbUID)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testMultiSaversXDMF(self, dbUID='GroupId5'):
        self._testMultiSaversDB(XDMFHandler, self.hdf5Args, self.hdf5Kwargs, dbUID=dbUID)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testXDMFSaving(self, dbUID='GroupId6'):
        self._testDBSaving(HDF5Handler, 'HDF5')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testXDMFGridSaving(self, dbUID='GroupID7'):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testXDMFGridSaving', suffix='')
        self.createdFiles.add(dbPath)

        rs = ResultSelector(self.newSim)
        saver1 = rs.ALL(Compartment).ALL(Species).Count
        saver1 <<= rs.ALL(Patch).ALL(Species).Count
        saver1 <<= rs.comp1.CC.sus2.Count

        tets = self.newGeom.tets
        tris1 = self.newGeom.patch.tris
        tris2 = self.newGeom.patch2.tris
        saver2 = rs.TETS(tets).S1.Count << rs.TRIS(tris1 + tris2).ExS1S2.Count
        if self.useEField:
            saver2 = saver2 << rs.VERTS(tris1.verts).V

        saver3 = rs.TETS(tets[0:len(tets)//2]).S2.Count << rs.TRIS(tris1).Ex.Count << rs.diffb.S1.DiffusionActive

        if self.useDist:
            saver4 = rs.ALL(Compartment).ALL(Species).Count
        else:
            saver4 = rs.ALL(ROI).ALL(Species).Count

        if self.useEField:
            saver5 = rs.TRIS(tris1).Chan1_Ohm_I.I
        else:
            saver5 = rs.TRIS(tris1).ExS1S2.Count

        if self.useDist:
            saver6 = rs.TRIS(tris1).Ex.Count
        else:
            saver6 = rs.TRIS(tris1).Area

        # Combined values accross geometrical elements
        saver7 = rs.SUM(rs.TRIS(tris1).Ex.Count)

        # Combined values on same geometrical elements
        saver8 = rs.TETS().S1.Count + rs.TETS().S2.Count

        self.newSim.toSave(saver1, saver2, saver3, saver4, saver5, saver6, saver7, saver8, dt=self.deltaT)

        with XDMFHandler(dbPath) as dbh:
            self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            for rid in range(self.nbRuns):
                self.newSim.newRun()
                self.init_API2_sim(self.newSim)
                self.newSim.run(self.shortEndTime)

            self.newSim.toDB(dbh, f'{dbUID}_2', val1=10, val2=20)
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

        with XDMFHandler(dbPath) as dbh:
            self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testXDMFColumnRemap(self, dbUID='GroupID8'):
        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testXDMFColumnRemap', suffix='')
        self.createdFiles.add(dbPath)
        if self.useDist:
            # Need to synchronize the file path prefix across ranks
            import mpi4py
            dbPath = mpi4py.MPI.COMM_WORLD.bcast(dbPath, root=0)

        S1, S2 = self.newSim.model.S1, self.newSim.model.S2

        rs = ResultSelector(self.newSim)

        saver1 = rs.TETS().LIST(S1, S2).Count
        saver2 = rs.TETS().S1.Count + rs.TETS().S2.Count
        saver3 = rs.SUM(rs.TETS().S1.Count) << rs.SUM(2 * rs.TETS().S2.Count)

        self.newSim.toSave(saver1, saver2, saver3, dt=self.deltaT)

        vals1 = list(range(len(self.newSim.geom.tets)))
        vals2 = list(reversed(vals1))
        with XDMFHandler(dbPath) as dbh:
            self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.TETS().S1.Count = vals1
            self.newSim.TETS().S2.Count = vals2
            self.newSim.run(self.shortEndTime)
            self.createdFiles |= set(dbh._getFilePaths())

        if self.useDist:
            # Need an MPI barrier to make sure the files are closed in all ranks
            mpi4py.MPI.COMM_WORLD.Barrier()

        # Load data
        if MPI._shouldWrite:
            with XDMFHandler(dbPath) as dbh:
                saver1, saver2, saver3 = dbh.get().results
                self.assertEqual(list(saver1.data[0,0,0::2]), vals1)
                self.assertEqual(list(saver1.data[0,0,1::2]), vals2)
                self.assertEqual(list(saver2.data[0,0,:]), list(np.array(vals1) + np.array(vals2)))
                self.assertEqual(list(saver3.data[0,0,:]), [sum(vals1), 2 * sum(vals2)])

        if self.useDist:
            # Need an MPI barrier to prevent early removal of files
            mpi4py.MPI.COMM_WORLD.Barrier()

    def testSimPathFunctionSetting(self):
        sim = self.newSim
        self.init_API2_sim(sim)
        self.specInds = {spec: i for i, spec in enumerate(sim.model.ALL(Species))}

        sim.TETS().ALL(Species).Count = lambda tet, spec: self.specInds[spec] + len(self.specInds) * tet.idx

        for tet in sim.geom.tets:
            for spec in tet.ALL(Species):
                self.assertEqual(sim.TET(tet).LIST(spec).Count, self.specInds[spec] + len(self.specInds) * tet.idx)

    @unittest.skip('Not needed here')
    def testXDMFWithoutMesh(self):
        pass

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testXDMFWithoutMPI(self, dbUID='GroupID'):
        """Check that XDMFHandler does not trigger any MPI initialization or any mpi4py import"""

        _, dbPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testXDMFWithoutMPISaving', suffix='')
        self.createdFiles.add(dbPath)

        # Setup traps
        steps.stepslib.mpiInit = lambda: self.assertTrue(False)
        sys.modules['mpi4py.MPI'] = None

        rs = ResultSelector(self.newSim)
        saver = rs.TETS().LIST('S1', 'S2').Count
        self.newSim.toSave(saver, dt=self.deltaT)

        with XDMFHandler(dbPath) as dbh:
            self.newSim.toDB(dbh, dbUID, val1=1, val2=2)
            self.createdFiles |= set(dbh._getFilePaths())
            self.newSim.newRun()
            self.init_API2_sim(self.newSim)
            self.newSim.ALL(Compartment).LIST('S1', 'S2').Count = 100
            self.newSim.run(self.shortEndTime)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(ResultSelectorTests))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(SimDataSaving))
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetSimDataSaving))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
