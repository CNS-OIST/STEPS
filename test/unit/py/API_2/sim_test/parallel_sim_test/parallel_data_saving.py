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

""" Unit tests for parallel data saving."""

import unittest
import numpy

from steps import interface

from steps.geom import *
from steps.rng import *
from steps.sim import *

import steps.API_1.mpi.solver as ssolver
import steps.API_1.mpi as smpi
import steps.API_1.rng as srng
import steps.API_1.utilities.geom_decompose as gd

from ..test_dataSaving import TetSimDataSaving
from .. import base_model

class TetOpSplitSimDataSaving(TetSimDataSaving):
    """Test data access, setting, and saving with tetmeshes."""

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        part = LinearMeshPartition(ngeom, 1, MPI.nhosts, 1)
        nsim = Simulation('TetOpSplit', nmdl, ngeom, nrng, MPI.EF_DEFAULT, part)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        orng=srng.create('mt19937',512)
        orng.initialize(self.seed)
        tet_hosts = gd.linearPartition(ogeom, [1, smpi.nhosts, 1])
        tri_hosts = gd.partitionTris(ogeom, tet_hosts, (self.ptris + self.ptris2).indices)
        if smpi.rank == 0:
            gd.validatePartition(ogeom, tet_hosts)
        osim = ssolver.TetOpSplit(omdl, ogeom, orng, ssolver.EF_DEFAULT, tet_hosts, tri_hosts)
        osim.setEfieldDT(self.efielddt)
        return osim


class VesRaftRDEFSimDataSaving(base_model.VesTestModelFramework, TetSimDataSaving):
    """Test data access, setting, and saving with vesicle and raft simulations."""

    def setUp(self):
        base_model.VesTestModelFramework.setUp(self)
        TetSimDataSaving.setUp(self, False)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('TetVesicle', nmdl, ngeom, nrng, False)

        self.one_time_init_API2_sim(nsim)

        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        orng=srng.create('mt19937',512)
        orng.initialize(self.seed)
        osim = ssolver.TetVesicle(omdl, ogeom, orng, False)
        osim.setOutputSync(True, 0)

        self.one_time_init_API1_sim(osim)

        return osim

    def testAddAndDeleteVesicle(self):
        sim = self.newSim
        mdl = self.newMdl

        self.assertEqual(sim.comp1.ves1.Count, 0)

        ves1_1 = sim.comp1.addVesicle(mdl.ves1)

        self.assertEqual(sim.VESICLE(ves1_1).Compartment, 'comp1')
        self.assertEqual(ves1_1.Compartment, 'comp1')
        self.assertTrue(ves1_1.exists())
        self.assertEqual(sim.comp1.ves1.Count, 1)

        ves1_2 = sim.comp1.addVesicle(mdl.ves1)

        self.assertEqual(sim.VESICLE(ves1_2).Compartment, 'comp1')
        self.assertEqual(ves1_2.Compartment, 'comp1')
        self.assertTrue(ves1_2.exists())
        self.assertEqual(sim.comp1.ves1.Count, 2)

        ves1_1.delete()

        self.assertFalse(ves1_1.exists())
        self.assertEqual(sim.comp1.ves1.Count, 1)

        ves1_2.delete()

        self.assertFalse(ves1_2.exists())
        self.assertEqual(sim.comp1.ves1.Count, 0)

        sim.comp1.ves1.Count = 5

        self.assertEqual(sim.comp1.ves1.Count, 5)

    def testAddAndDeleteRaft(self):
        sim = self.newSim
        mdl = self.newMdl
        tri = sim.patch.tris[0]

        self.assertEqual(sim.TRI(tri).raft1.Count, 0)
        self.assertEqual(sim.patch.raft1.Count, 0)

        raft1_1 = sim.TRI(tri).addRaft(mdl.raft1)

        self.assertTrue(raft1_1.exists())
        self.assertEqual(sim.RAFT(raft1_1).Patch, 'patch')
        self.assertEqual(raft1_1.Patch, 'patch')
        self.assertEqual(sim.TRI(tri).raft1.Count, 1)
        self.assertEqual(sim.patch.raft1.Count, 1)

        # raft1_1.delete()
        # self.assertFalse(raft1_1.exists())

        sim.patch.raft1.Count = 2

        self.assertEqual(sim.patch.raft1.Count, 2)

        sim.patch.raft1.Count = 0
        self.assertEqual(sim.patch.raft1.Count, 0)

        sim.TRI(tri).raft1.Count = 3

        self.assertEqual(sim.TRI(tri).raft1.Count, 3)
        self.assertEqual(sim.patch.raft1.Count, 3)

    def testMoveVesicles(self):
        sim = self.newSim
        mdl = self.newMdl

        ves1_1 = sim.comp1.addVesicle(mdl.ves1)
        ves1_2 = sim.comp1.addVesicle(mdl.ves1)

        pos1 = ves1_1.Pos
        pos2 = ves1_2.Pos

        self.assertNotEqual(tuple(pos1), tuple(pos2))

        # TODO change setCompSingleVesiclePos -> setSingleVesiclePos
        # ves1_1.Pos = ves1_2.Pos

        # self.assertEqual(tuple(ves1_1.Pos), tuple(pos1))
        # self.assertEqual(tuple(ves1_2.Pos), tuple(pos2))

        ves1_1.setPos(ves1_2.Pos, force=True)

        self.assertEqual(tuple(ves1_1.Pos), tuple(pos2))
        self.assertEqual(tuple(ves1_2.Pos), tuple(pos1))

    def testMoveRafts(self):
        sim = self.newSim
        mdl = self.newMdl
        tri1 = sim.patch.tris[0]
        tri2 = sim.patch.tris[1]

        raft1_1 = sim.TRI(tri1).addRaft(mdl.raft1)
        raft1_2 = sim.TRI(tri2).addRaft(mdl.raft1)

        pos1 = raft1_1.Pos
        pos2 = raft1_2.Pos

        self.assertNotEqual(tuple(pos1), tuple(pos2))

        # TODO Setting raft position not implemented yet

    def testVesicleLists(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.comp1.ves1.Count = 5

        vesicles = VesicleList(sim.comp1.ves1)

        positions = vesicles.Pos

        self.assertEqual(len(positions), 5)
        self.assertEqual(len(set(map(tuple, positions))), 5)

        self.assertEqual(vesicles.indices, sim.comp1.ves1.Indices)

    def testRaftLists(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.patch.raft1.Count = 5

        rafts = RaftList(sim.patch.raft1)

        positions = rafts.Pos

        self.assertEqual(len(positions), 5)
        # Rafts can overlap
        # self.assertEqual(len(set(map(tuple, positions))), 5)

        self.assertEqual(rafts.indices, sim.patch.raft1.Indices)

    def testLinkSpecLists(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.newRun()

        sim.comp1.ves1.Count = 10

        veslst = VesicleList(sim.comp1.ves1)

        sim.VESICLES(veslst('surf')).S1.Count = 10

        sim.run(1)

        ves2lsDct = {}
        for ves in veslst:
            L1lst = LinkSpecList(sim.VESICLE(ves)('surf').L1)
            allPos = [ls.Pos for ls in L1lst]

            self.assertEqual(L1lst.indices, sim.VESICLE(ves)('surf').L1.Indices)
            self.assertEqual(allPos, sim.VESICLE(ves)('surf').L1.Pos)

            posDct = sim.VESICLE(ves)('surf').LINKSPECS(mdl.L1).Pos
            self.assertEqual(list(posDct.values()), allPos)
            ves2lsDct[ves.idx] = posDct

        autoDct = sim.comp1.VESICLES(mdl.ves1)('surf').LINKSPECS(mdl.L1).Pos
        self.assertEqual(ves2lsDct, autoDct)

    def testSetSingleRaftCounts(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.patch.raft1.Count = 2

        raftlst = RaftList(sim.patch.raft1)
        raft1_1, raft1_2 = raftlst

        self.assertEqual(sim.patch.raft1.S1.Count, 0)
        self.assertEqual(sim.RAFT(raft1_1).S1.Count, 0)
        self.assertEqual(sim.RAFT(raft1_2).S1.Count, 0)
        self.assertEqual(sim.RAFTS(raftlst).S1.Count, [0, 0])

        raft1_1.S1.Count = 10

        self.assertEqual(sim.patch.raft1.S1.Count, 10)
        self.assertEqual(sim.RAFT(raft1_1).S1.Count, 10)
        self.assertEqual(sim.RAFT(raft1_2).S1.Count, 0)
        self.assertEqual(sim.RAFTS(raftlst).S1.Count, [10, 0])

        raft1_2.S1.Count = 15

        self.assertEqual(sim.patch.raft1.S1.Count, 25)
        self.assertEqual(sim.RAFT(raft1_1).S1.Count, 10)
        self.assertEqual(sim.RAFT(raft1_2).S1.Count, 15)
        self.assertEqual(sim.RAFTS(raftlst).S1.Count, [10, 15])

        sim.RAFTS(raftlst).S2.Count = 20

        self.assertEqual(sim.patch.raft1.S2.Count, 40)
        self.assertEqual(sim.RAFT(raft1_1).S2.Count, 20)
        self.assertEqual(sim.RAFT(raft1_2).S2.Count, 20)
        self.assertEqual(sim.RAFTS(raftlst).S2.Count, [20, 20])

        sim.RAFTS(raftlst).S2.Count = [31, 35]

        self.assertEqual(sim.patch.raft1.S2.Count, 66)
        self.assertEqual(sim.RAFT(raft1_1).S2.Count, 31)
        self.assertEqual(sim.RAFT(raft1_2).S2.Count, 35)
        self.assertEqual(sim.RAFTS(raftlst).S2.Count, [31, 35])

    def testSetSingleVesicleCounts(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.comp1.ves1.Count = 2

        # Test both surface and inside
        for loc in ['surf', 'in']:
            veslst = VesicleList(sim.comp1.ves1)
            ves1_1, ves1_2 = veslst

            self.assertEqual(sim.comp1.ves1(loc).S1.Count, 0)
            self.assertEqual(sim.VESICLE(ves1_1)(loc).S1.Count, 0)
            self.assertEqual(sim.VESICLE(ves1_2)(loc).S1.Count, 0)
            self.assertEqual(sim.VESICLES(veslst)(loc).S1.Count, [0, 0])

            ves1_1(loc).S1.Count = 10

            self.assertEqual(sim.comp1.ves1(loc).S1.Count, 10)
            self.assertEqual(sim.VESICLE(ves1_1)(loc).S1.Count, 10)
            self.assertEqual(sim.VESICLE(ves1_2)(loc).S1.Count, 0)
            self.assertEqual(sim.VESICLES(veslst)(loc).S1.Count, [10, 0])

            ves1_2(loc).S1.Count = 15

            self.assertEqual(sim.comp1.ves1(loc).S1.Count, 25)
            self.assertEqual(sim.VESICLE(ves1_1)(loc).S1.Count, 10)
            self.assertEqual(sim.VESICLE(ves1_2)(loc).S1.Count, 15)
            self.assertEqual(sim.VESICLES(veslst)(loc).S1.Count, [10, 15])

            sim.VESICLES(veslst)(loc).S2.Count = 20

            self.assertEqual(sim.comp1.ves1(loc).S2.Count, 40)
            self.assertEqual(sim.VESICLE(ves1_1)(loc).S2.Count, 20)
            self.assertEqual(sim.VESICLE(ves1_2)(loc).S2.Count, 20)
            self.assertEqual(sim.VESICLES(veslst)(loc).S2.Count, [20, 20])

            sim.VESICLES(veslst)(loc).S2.Count = [31, 35]

            self.assertEqual(sim.comp1.ves1(loc).S2.Count, 66)
            self.assertEqual(sim.VESICLE(ves1_1)(loc).S2.Count, 31)
            self.assertEqual(sim.VESICLE(ves1_2)(loc).S2.Count, 35)
            self.assertEqual(sim.VESICLES(veslst)(loc).S2.Count, [31, 35])

    def testGetSetVesicleSpecPos(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.comp1.ves1.Count = 2

        veslst = VesicleList(sim.comp1.ves1)
        ves1_1, ves1_2 = veslst

        ves1_1('surf').S1.Count = 10

        positions = ves1_1('surf').S1.Pos
        center = numpy.array(ves1_1.Pos)

        self.assertEqual(len(positions), 10)

        for pos in positions:
            self.assertEqual(len(pos), 3)
            self.assertAlmostEqual(numpy.linalg.norm(numpy.array(pos) - center), mdl.ves1.Diameter / 2, places=15)

        # set spherical spec pos

        ves1_1('surf').S1.PosSpherical = [0, 0]

        positions = [tuple(pos) for pos in ves1_1('surf').S1.Pos]

        self.assertEqual(len(set(positions)), 1)
        x, y, z = positions[0]
        cx, cy, cz = center

        self.assertAlmostEqual(x, cx, places=15)
        self.assertAlmostEqual(y, cy, places=15)
        self.assertAlmostEqual(z, cz + mdl.ves1.Diameter / 2, places=15)

        ves1_1('surf').S1.PosSpherical = [numpy.pi / 2, numpy.pi / 2]
        positions = [tuple(pos) for pos in ves1_1('surf').S1.Pos]
        x, y, z = positions[0]
        self.assertAlmostEqual(x, cx, places=15)
        self.assertAlmostEqual(y, cy + mdl.ves1.Diameter / 2, places=15)
        self.assertAlmostEqual(z, cz, places=15)

    def testGetSingleVesOverlapTets(self):
        sim = self.newSim
        mdl = self.newMdl
        mesh = self.newGeom

        sim.comp1.ves1.Count = 2

        veslst = VesicleList(sim.comp1.ves1)
        ves1_1, ves1_2 = veslst

        tetlst = TetList(ves1_1.OverlapTets, mesh=mesh)
        x, y, z = ves1_1.Pos

        self.assertGreaterEqual(len(tetlst), 1)

        self.assertTrue(mesh.tets[x, y, z] in tetlst)

    def testsetTetVesicleDcst(self):
        sim = self.newSim
        mdl = self.newMdl
        mesh = self.newGeom

        sim.TET(mesh.comp1.tets[0]).ves1.Dcst = mdl.ves1.Dcst / 2

    def testsetLinkSPecSDiffD(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.ves1('surf').L1.SDiffD = 1e-5

    def testGetSingleRaftImmobility(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.patch.raft1.Count = 2

        raftlst = RaftList(sim.patch.raft1)
        raft1_1, raft1_2 = raftlst

        self.assertEqual(raft1_1.Immobility, 0)
        self.assertEqual(raft1_2.Immobility, 0)

        self.assertEqual(raftlst.Immobility, [0, 0])

    def testGetSingleVesicleImmobility(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.comp1.ves1.Count = 2

        veslst = VesicleList(sim.comp1.ves1)
        ves1_1, ves1_2 = veslst

        self.assertEqual(ves1_1.Immobility, 0)
        self.assertEqual(ves1_2.Immobility, 0)

        self.assertEqual(veslst.Immobility, [0, 0])

    def testLinkSpecCountAndPos(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.newRun()

        sim.comp1.ves1.Count = 10

        veslst = VesicleList(sim.comp1.ves1)

        sim.VESICLES(veslst('surf')).S1.Count = 10
        sim.ves1sreac1.K = 0, 0

        sim.run(1)

        counts = sim.VESICLES(veslst)('surf').S1.Count

        self.assertTrue(any(cnt < 10 for cnt in counts))

        counts2 = sim.VESICLES(veslst)('surf').L1.Count

        self.assertTrue(all(counts2[idx] > 0 if cnt < 10 else counts2[idx] == 0 for idx, cnt in enumerate(counts)))

        allPositions = sim.VESICLES(veslst)('surf').L1.Pos

        self.assertGreater(len(allPositions), 0)

        for positions, ves in zip(allPositions, veslst):
            center = numpy.array(ves.Pos)
            for pos in positions:
                self.assertAlmostEqual(numpy.linalg.norm(numpy.array(pos) - center), mdl.ves1.Diameter / 2, places=15)

    def testRaftEndocytosisKcstChanges(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.newRun()

        sim.patch.raft1.Count = 5

        for raft in RaftList(sim.patch.raft1):
            self.assertEqual(raft.rendo1.K, self.rendo1Kcst)
            raft.rendo1.K = self.rendo1Kcst * 10
            self.assertEqual(raft.rendo1.K, self.rendo1Kcst * 10)

        sim.rendo1.K = self.rendo1Kcst

        for raft in RaftList(sim.patch.raft1):
            self.assertEqual(raft.rendo1.K, self.rendo1Kcst)

    def testEndocytosisKcstChanges(self):
        sim = self.newSim

        sim.newRun()

        sim.patch.endoZone1.endo1.K = 2 * self.endo1Kcst

    # TODO Move this test to validation, it takes too long to be part of unit tests
    def _testUnbindingRates(self):
        sim = self.newSim
        mdl = self.newMdl

        sim.newRun()

        sim.comp1.ves2.Count = 20
        veslst = VesicleList(sim.comp1.ves2)

        self.assertEqual(len(veslst), 20)

        sim.VESICLES(veslst)('surf').S1.Count = 1000

        time = 10
        timeFrac = 0.1

        # The binding reaction is much faster than the unbinding so we assume it happened
        # completely at t1 = time * timeFrac. Since the binding reaction immobilizes the vesicles,
        # the number of formed links only goes down for the rest of the simulation, we can thus
        # estimate the rate at which it goes down from the number of broken links.

        sim.run(timeFrac * time)

        nbL1t1 = sum(sim.comp1.VESICLES(mdl.ves2)('surf').L1.Count.values()) / 2
        nbL2t1 = sum(sim.comp1.VESICLES(mdl.ves2)('surf').L2.Count.values())

        sim.run(time)

        nbL1t2 = sum(sim.comp1.VESICLES(mdl.ves2)('surf').L1.Count.values()) / 2
        nbL2t2 = sum(sim.comp1.VESICLES(mdl.ves2)('surf').L2.Count.values())

        # Estimate the unbinding reaction rate
        rateL1 = - numpy.log(nbL1t2 / nbL1t1) / (time * (1 - timeFrac))
        rateL2 = - numpy.log(nbL2t2 / nbL2t1) / (time * (1 - timeFrac))

        actualRateL1 = self.ves2unbind1Kcst + self.ves2unbind2Kcst
        actualRateL2 = self.ves2unbind3Kcst

        # Check that the estimated rates are within 25% of the actual rates
        self.assertLess(abs(rateL1 - actualRateL1) / actualRateL1, 0.25)
        self.assertLess(abs(rateL2 - actualRateL2) / actualRateL2, 0.25)

        # Check that the ratio of S2 to S3 matches the ratio of the rates within 25%
        nbS2 = sum(sim.comp1.VESICLES(mdl.ves2)('surf').S2.Count.values())
        nbS3 = sum(sim.comp1.VESICLES(mdl.ves2)('surf').S3.Count.values())

        ratio = self.ves2unbind1Kcst / self.ves2unbind2Kcst

        self.assertLess(abs(nbS2 / nbS3 - ratio) / ratio, 0.25)

    def _testDeltaTSavingParams(self):
        timepoints = numpy.arange(0, self.endTime + self.deltaT, self.deltaT)
        def oldSave(sim):
            res = []
            # TODO Add data to be saved
            res.append(sim.getCompVesicleSurfaceSpecCount('comp1', 'ves1', 'S1'))
            return res

        def newSave(rs):
            saver = rs.comp1.ves1('surf').S1.Count
            self.newSim.toSave(saver, dt=self.deltaT)
            return [saver]

        def newGetData(savers):
            saver = savers[0]
            saver.time[0]
            numpy.array(saver.time)
            saver.data[0,:]
            numpy.array(saver.data)
            return saver.data[:]

        def newRun(sim, rs):
            sim.run(self.endTime)

        explbls = [[None]]

        return (
            (timepoints, oldSave),
            (newSave, newGetData, newRun, explbls),
            [self.countRefVal]*2
        )
    

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TetOpSplitSimDataSaving, "test"))
    all_tests.append(unittest.makeSuite(VesRaftRDEFSimDataSaving, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
