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

import os
import tempfile
import unittest
import mpi4py.MPI

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

from . import base_model

def getUniqueTempPrefix(*args, **kwargs):
    _, pathPrefix = tempfile.mkstemp(*args, **kwargs)
    os.remove(pathPrefix)
    if MPI._nhosts > 1:
        # Need to synchronize across ranks
        pathPrefix = mpi4py.MPI.COMM_WORLD.bcast(pathPrefix, root=0)
    return pathPrefix

class WmCheckpoints(base_model.TestModelFramework):
    """Test checkpoint and restore methods for Wmdirect solver"""
    def setUp(self, callParent=True):
        if callParent:
            super().setUp()

        # Use true complexes for checkpoint tests
        self.newMdl = self.get_API2_Mdl(sas=self.useMesh)
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

    def _getCPFiles(self, pathPrefix):
        tmpDir, prefix = os.path.split(pathPrefix)
        cpfiles = []
        for name in os.listdir(tmpDir):
            if name.startswith(prefix):
                fullPath = os.path.join(tmpDir, name)
                if MPI._nhosts > 1:
                    *_, rank = name.split('_')
                    rank = int(rank)
                    self.assertLess(rank, MPI._nhosts)
                    self.assertGreaterEqual(rank, 0)
                    if rank != MPI._rank:
                        continue
                    noRankPath = os.path.join(tmpDir, name[:-len(f'_{rank}')])
                else:
                    noRankPath = fullPath
                cpfiles.append((noRankPath, fullPath))
        return cpfiles

    def _getAutoCPFiles(self, pathPrefix):
        tmpDir, prefix = os.path.split(pathPrefix)
        cpfiles = []
        for name in os.listdir(tmpDir):
            if name.startswith(prefix):
                fullPath = os.path.join(tmpDir, name)
                run, time, solver, *cp = name[(len(prefix)+1):].split('_')
                run = int(run)
                time = float(time)
                if MPI._nhosts > 1:
                    cp, rank = cp
                    rank = int(rank)
                    self.assertLess(rank, MPI._nhosts)
                    self.assertGreaterEqual(rank, 0)
                    if rank != MPI._rank:
                        continue
                    noRankPath = os.path.join(tmpDir, name[:-len(f'_{rank}')])
                else:
                    cp = cp[0]
                    noRankPath = fullPath
                self.assertEqual(cp, 'cp')
                cpfiles.append((run, time, noRankPath, fullPath))
        cpfiles.sort()
        return cpfiles

    def _cleanCPFiles(self, cpFiles):
        for *_, fullPath in cpFiles:
            os.remove(fullPath)

    def testBasicCheckpoints(self):
        cpPath = getUniqueTempPrefix(prefix=f'{self.__class__.__name__}basic')
        _, cpEmptyPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}empty')

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

        self._cleanCPFiles(self._getCPFiles(cpPath))

    def testAutoCheckpoints(self):
        pathPrefix = getUniqueTempPrefix(prefix=f'{self.__class__.__name__}auto')

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

        cpfiles = self._getAutoCPFiles(pathPrefix)

        self.assertEqual(len(cpfiles), self.nbRuns * (self.endTime / self.cpTime + 1))
        for r in range(self.nbRuns):
            self.assertEqual(sum(1 for run, time, *_ in cpfiles if run == r), self.endTime / self.cpTime + 1)

        for run, time, cpPath, cpFullPath in cpfiles:
            self.newSim.newRun()

            self.newSim.restore(cpPath)

            self.assertAlmostEqual(self.newSim.Time, time)

        self._cleanCPFiles(cpfiles)

    def testOnlyLast(self):
        pathPrefix = getUniqueTempPrefix(prefix=f'{self.__class__.__name__}onlyLast')

        self.newSim.autoCheckpoint(self.cpTime, pathPrefix, onlyLast=True)

        for r in range(self.nbRuns):
            self.newSim.newRun()
            self.newSim.run(self.endTime)

        cpfiles = self._getAutoCPFiles(pathPrefix)
        
        self.assertEqual(len(cpfiles), 1)
        self.assertEqual(cpfiles[0][0], self.nbRuns - 1)

        cpPath = cpfiles[0][2]

        self.newSim.newRun()

        self.newSim.restore(cpPath)

        self.assertAlmostEqual(self.newSim.Time, cpfiles[0][1])

        self._cleanCPFiles(cpfiles)
        
    def testStopAutoCheckpointing(self):
        pathPrefix = getUniqueTempPrefix(prefix=f'{self.__class__.__name__}stopAutoChkpt')

        for r in range(self.nbRuns):
            self.newSim.newRun()
            self.newSim.autoCheckpoint(self.cpTime, pathPrefix)
            self.newSim.run(self.endTime / 2)
            self.newSim.autoCheckpoint(None)
            self.newSim.run(self.endTime)

        cpfiles = self._getAutoCPFiles(pathPrefix)

        self.assertEqual(len(cpfiles), self.nbRuns * (self.endTime / 2 / self.cpTime + 1))
        for r in range(self.nbRuns):
            self.assertEqual(sum(1 for run, time, *_ in cpfiles if run == r), self.endTime / 2 / self.cpTime + 1)
            self.assertTrue(all(time <= self.endTime / 2 for run, time, *_ in cpfiles if run == r))

        for run, time, cpPath, cpFullPath in cpfiles:
            self.newSim.newRun()

            self.newSim.restore(cpPath)

            self.assertAlmostEqual(self.newSim.Time, time)

        self._cleanCPFiles(cpfiles)

    def testContinueRestoredRun(self):
        pathPrefix = getUniqueTempPrefix(prefix=f'{self.__class__.__name__}ContinueRestoredRun')

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

        self._cleanCPFiles(self._getCPFiles(pathPrefix))

    def testReproducibleCheckpointing(self, plot=False):
        pathPrefix = getUniqueTempPrefix(prefix=f'{self.__class__.__name__}ReproducibleChekpointing')

        def getSelectors(sim):
            rs = ResultSelector(sim)
            selectors = []
            # Counts
            selectors.append(rs.ALL(Compartment, Patch).ALL(Species).Count)
            # Extents
            selectors.append(rs.ALL(Compartment).ALL(Reaction).Extent)
            # Membrane potential
            if self.useEField and not self.useDist:
                verts = self.newGeom.membrane.tris.verts
                selectors.append(rs.SUM(rs.VERTS(verts).V) / len(verts))

            return selectors

        selectors = getSelectors(self.newSim)
        self.newSim.toSave(*selectors, dt=self.deltaT)

        self.newSim.newRun()

        self.init_API2_sim(self.newSim)
        self.newSim.run(self.endTime / 2)
        self.newSim.checkpoint(pathPrefix)
        self.newSim.run(self.endTime)

        # Continue run in same simulation
        self.newSim.newRun()
        self.newSim.restore(pathPrefix)
        self.newSim.run(self.endTime)

        # Continue run in different simulation
        self.newSim = self._get_API2_Sim(self.newMdl, self.newGeom)
        selectors2 = getSelectors(self.newSim)
        self.newSim.toSave(*selectors2, dt=self.deltaT)

        self.newSim.newRun()
        self.newSim.restore(pathPrefix)
        self.newSim.run(self.endTime)

        if MPI._shouldWrite:
            for sel1, sel2 in zip(selectors, selectors2):
                n = len(sel1.data[0]) - len(sel1.data[1])
                reference = sel1.data[0, n:, :]
                restored1 = sel1.data[1, :, :]
                restored2 = sel2.data[0, :, :]

                # Plotting code in case manual inspection is needed
                if plot:
                    from matplotlib import pyplot as plt
                    for i in range(reference.shape[-1]):
                        if (reference[:, i] != restored1[:, i]).any() or (reference[:, i] != restored2[:, i]).any():
                            plt.plot(sel1.time[1], reference[:, i], label='Reference run')
                            plt.plot(sel1.time[1], restored1[:, i], '--', label='Restored run, same simulation')
                            plt.plot(sel1.time[1], restored2[:, i], '--', label='Restored run, new simulation')
                            plt.legend()
                            plt.title(sel1.labels[i])
                            plt.show()

                # Check that all runs have exactly the same traces
                self.assertTrue((reference == restored1).all())
                self.assertTrue((reference == restored2).all())

        self._cleanCPFiles(self._getCPFiles(pathPrefix))


class TetCheckpoints(base_model.TetTestModelFramework, WmCheckpoints):
    """Test checkpoint and restore methods for Tetexact solver"""

    def setUp(self):
        base_model.TetTestModelFramework.setUp(self)
        self.useEField = True
        self.useVesicle = False

        WmCheckpoints.setUp(self, False)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Tetexact', nmdl, ngeom, nrng, True)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _setInitValsToCurrentOnes(self):
        WmCheckpoints._setInitValsToCurrentOnes(self)

        sim = self.newSim

        if self.useEField:
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

