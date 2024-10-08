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

""" Unit tests for distmesh class and related methods."""

import functools
import numpy
import operator
import os
import tempfile
import unittest

from steps import interface

from steps.geom import *
from steps.rng import *
from steps.sim import *

from . import test_tetMesh

FILEDIR = os.path.dirname(os.path.abspath(__file__))


class distTetMeshTests(test_tetMesh.tetMeshTests):
    def setUpMeshes(self):
        # Deactivate tests on features that are not yet implemented
        self._distMesh = True
        self.refMesh = TetMesh.LoadAbaqus(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets.inp'), 1e-6)

        self.mesh = DistMesh(os.path.join(FILEDIR, '../../../../mesh/box.msh'))
        self.mesh2 = DistMesh(os.path.join(FILEDIR, '../../../../mesh/box.msh'))

        self.mesh3 = DistMesh(os.path.join(FILEDIR, 'meshes', 'tagged.msh'))

        self.elemListNames = ['tets', 'tris', 'bars', 'verts']

    def getSimulation(self):
        rng = RNG('mt19937', 512, 7233)
        return Simulation('DistTetOpSplit', self.model, self.mesh, rng)

    def gather(self, lst, func=lambda a, b: a + b):
        lists = self.mesh3._comm.gather(lst, root=0)
        if lists is not None:
            return functools.reduce(func, lists)
        else:
            return None

    def allGather(self, lst, func=lambda a, b: a + b):
        lists = self.mesh3._comm.allgather(lst)
        return functools.reduce(func, lists)

    def testGroupBlocksLoading(self):
        mesh = self.mesh3

        self.assertSetEqual(set(mesh.vertGroups.keys()), set(['NODEGROUP1', 'NODEGROUP2', 'NODEGROUP3']))
        self.assertSetEqual(set(mesh.triGroups.keys()), set(['TRIGROUP1', 'TRIGROUP2', 'TRIGROUP3']))
        self.assertSetEqual(set(mesh.tetGroups.keys()), set(['TETGROUP1', 'TETGROUP2', 'TETGROUP3']))

        tets = set(mesh.tetGroups['TETGROUP1'] + mesh.tetGroups['TETGROUP2'] + mesh.tetGroups['TETGROUP3'])
        self.assertSetEqual(tets, set(mesh.tets))

        self.assertEqual(set(mesh.triGroups['TRIGROUP1']), set(mesh.tetGroups['TETGROUP1'].surface))
        self.assertEqual(set(mesh.triGroups['TRIGROUP2']), set(mesh.tetGroups['TETGROUP2'].surface))
        self.assertEqual(set(mesh.triGroups['TRIGROUP3']), set(mesh.surface - mesh.tetGroups['TETGROUP1'].surface))

        self.assertEqual(len(set(mesh.vertGroups['NODEGROUP1'] & mesh.vertGroups['NODEGROUP3'])), 1)
        self.assertEqual(len(set(mesh.vertGroups['NODEGROUP1'] & mesh.vertGroups['NODEGROUP2'])), 0)

        with self.assertRaises(NotImplementedError):
            mesh.tetGroups['test'] = TetList([1, 2, 3], mesh=mesh)

        with self.assertRaises(NotImplementedError):
            del mesh.tetGroups['TETGROUP1']
        tets = set(mesh.tetGroups['TETGROUP1'] + mesh.tetGroups['TETGROUP2'] + mesh.tetGroups['TETGROUP3'])
        self.assertSetEqual(tets, set(mesh.tets))

    def testLocalLists(self):
        mesh = self.mesh3

        for useAsLocal in [False, True]:
            for attr in self.elemListNames:
                glob = getattr(mesh, attr)
                if useAsLocal:
                    with mesh.asLocal():
                        local = getattr(mesh, attr)
                else:
                    local = glob.toLocal()
                globLocal = local.toGlobal()

                self.assertTrue(local.isLocal())
                self.assertFalse(glob.isLocal())
                self.assertFalse(globLocal.isLocal())

                self.assertFalse(local == glob)
                self.assertFalse(local[0] == glob[0])

                self.assertEqual(len(local), len(globLocal))

                if MPI.nhosts > 1:
                    self.assertLess(len(local), len(glob))
                else:
                    self.assertEqual(len(local), len(glob))
                lengths = mesh._comm.gather(len(local), root=0)
                if lengths is not None:
                    self.assertEqual(sum(lengths), len(glob))

                indices = self.gather(globLocal.indices)
                if indices is not None:
                    self.assertEqual(set(indices), set(glob.indices))

    def testLocalGlobalContextManager(self):
        mesh = self.mesh3
        with mesh:
            tets = mesh.tets
            surf = mesh.surface
        with mesh.asLocal():
            localTets = mesh.tets
            localSurf = mesh.surface
            self.assertCountEqual(localTets.combineWithOperator(operator.or_), tets)
            self.assertCountEqual(localSurf.combineWithOperator(operator.or_), surf)
        with mesh.asGlobal():
            self.assertEqual(mesh.tets, tets)
            self.assertEqual(mesh.surface, surf)

        # No nesting
        with self.assertRaises(Exception):
            with mesh:
                with mesh.asLocal():
                    tets = mesh.tets

        with self.assertRaises(Exception):
            with mesh:
                with mesh.asGlobal():
                    tets = mesh.tets

    def testNonOwnedLocalLists(self):
        mesh = self.mesh3
        for attr in self.elemListNames:
            with mesh:
                lst = getattr(mesh, attr)
            with mesh.asLocal():
                ownedLst = getattr(mesh, attr)
            with mesh.asLocal(owned=False):
                ghostLst = getattr(mesh, attr)
            ownedLst2 = lst.toLocal()
            ghostLst2 = lst.toLocal(owned=False)
            self.assertCountEqual(ownedLst2, ownedLst)
            self.assertCountEqual(ghostLst2, ghostLst)
            self.assertCountEqual(ghostLst.toLocal(owned=True), ownedLst)
            self.assertCountEqual(ghostLst.toLocal(owned=False), ghostLst)
            self.assertLess(len(ownedLst), len(ghostLst))
            self.assertEqual(len(ownedLst - ghostLst), 0)
            self.assertGreater(len(ghostLst - ownedLst), 0)

            for locElem in ghostLst - ownedLst:
                self.assertIsNone(locElem.toLocal())
            for locElem in ownedLst:
                self.assertIsNotNone(locElem.toGlobal().toLocal())
            for locElem in ghostLst:
                self.assertIsNotNone(locElem.toGlobal().toLocal(owned=False))

            allLst1 = ownedLst.combineWithOperator(operator.or_)
            allLst2 = ghostLst.combineWithOperator(operator.or_)
            self.assertCountEqual(allLst1, lst)
            self.assertCountEqual(allLst2, lst)

    def testCombineWithOperator_n2(self):
        mesh = self.mesh3
        with mesh.asLocal(owned=False):
            lstAll = mesh.tets
        with mesh.asLocal():
            lstOwned = mesh.tets

        and_ = lstAll.combineWithOperator(operator.and_)
        ghosts = (lstAll - lstOwned).combineWithOperator(operator.or_)

        self.assertGreater(len(and_), 0)
        self.assertCountEqual(and_, ghosts)

        or_ = lstAll.combineWithOperator(operator.or_)

        self.assertCountEqual(or_, mesh.tets)

        addOwned = lstOwned.combineWithOperator(operator.add)
        addAll = lstAll.combineWithOperator(operator.add)

        self.assertCountEqual(addOwned, mesh.tets)
        self.assertGreater(len(addAll), len(mesh.tets))
        for g in ghosts:
            self.assertEqual(addAll.indices.count(g.idx), 2)

    def testLocalTetByPoint(self):
        mesh = self.mesh3

        points = [tet.center for tet in mesh.tets[::10]]
        for i, point in enumerate(points):
            tet = mesh.tets[point]
            self.assertFalse(tet.isLocal())
            self.assertEqual(tet.idx, i * 10)

            with mesh.asLocal():
                try:
                    localTet = mesh.tets[point]
                except KeyError:
                    localTet = None
                # Check that only one rank is not None
                isNone = mesh._comm.gather(localTet is None, root=0)
                if MPI.rank == 0:
                    self.assertEqual(isNone.count(False), 1)
                if localTet is not None:
                    self.assertTrue(localTet.isLocal())
                    self.assertEqual(localTet.idx, tet.toLocal().idx)
                    self.assertEqual(localTet.toGlobal().idx, tet.idx)
                    self.assertTrue(localTet.containsPoint(point))

    def testMeshProperties(self):
        mesh = self.mesh3

        # surface
        surf = mesh.surface
        with mesh.asLocal():
            localSurf = mesh.surface
            self.assertTrue(localSurf.isLocal())

            surf2 = self.gather(localSurf.toGlobal().indices)
            if surf2 is not None:
                self.assertEqual(set(surf2), set(surf.indices))
        self.assertEqual(set(surf), set(mesh.surface))

        # bbox
        bbox = mesh.bbox
        with mesh.asLocal():
            localMin = mesh.bbox.min
            localMax = mesh.bbox.max
            for l, g in zip(localMin, bbox.min):
                self.assertGreaterEqual(l, g)
            for l, g in zip(localMax, bbox.max):
                self.assertLessEqual(l, g)
            totalMin = self.gather(localMin, lambda a,b: [min(x,y) for x,y in zip(a,b)])
            totalMax = self.gather(localMax, lambda a,b: [max(x,y) for x,y in zip(a,b)])
            if totalMin is not None:
                self.assertEqual(list(totalMin), list(bbox.min))
            if totalMax is not None:
                self.assertEqual(list(totalMax), list(bbox.max))
        self.assertEqual(list(bbox.min), list(mesh.bbox.min))
        self.assertEqual(list(bbox.max), list(mesh.bbox.max))

        with mesh.asLocal(owned=False):
            with self.assertWarns(UserWarning):
                self.assertEqual(list(mesh.bbox.min), list(localMin))
                self.assertEqual(list(mesh.bbox.max), list(localMax))

        # Vol
        vol = mesh.Vol
        with mesh.asLocal():
            localVol = mesh.Vol
            self.assertLessEqual(localVol, vol)
            totVol = self.gather(localVol)
            if totVol is not None:
                self.assertEqual(totVol, vol)
        self.assertEqual(vol, mesh.Vol)

        with mesh.asLocal(owned=False):
            with self.assertWarns(UserWarning):
                self.assertEqual(mesh.Vol, localVol)

        # elem groups
        tets = mesh.tetGroups['TETGROUP3']
        tris = mesh.triGroups['TRIGROUP3']
        verts = mesh.vertGroups['NODEGROUP3']
        with mesh.asLocal():
            localTets = mesh.tetGroups['TETGROUP3']
            self.assertLessEqual(len(localTets), len(tets))
            totTets = self.gather(localTets.toGlobal().indices)
            if totTets is not None:
                self.assertEqual(set(totTets), set(tets.indices))

            localTris = mesh.triGroups['TRIGROUP3']
            self.assertLessEqual(len(localTris), len(tris))
            totTris = self.gather(localTris.toGlobal().indices)
            if totTris is not None:
                self.assertEqual(set(totTris), set(tris.indices))

            localVerts = mesh.vertGroups['NODEGROUP3']
            self.assertLessEqual(len(localVerts), len(verts))
            totVerts = self.gather(localVerts.toGlobal().indices)
            if totVerts is not None:
                self.assertEqual(set(totVerts), set(verts.indices))
        self.assertEqual(set(tets), set(mesh.tetGroups['TETGROUP3']))
        self.assertEqual(set(tris), set(mesh.triGroups['TRIGROUP3']))
        self.assertEqual(set(verts), set(mesh.vertGroups['NODEGROUP3']))

    def testTetListProperties(self):
        mesh = self.mesh3

        # surface
        meshSurf = mesh.surface
        self.assertEqual(set(meshSurf), set(mesh.tets.surface))

    def testDistReferences(self):
        mesh = self.mesh3

        for attr in self.elemListNames[1:]:
            elems = getattr(mesh, attr)
            for i, elem in enumerate(elems):
                self.assertFalse(elem.isLocal())
                localElem = elem.toLocal()
                self.assertEqual(elem.toGlobal(), elem)
                if localElem is not None:
                    localElem2 = localElem.toLocal()
                    self.assertTrue(localElem2.isLocal())
                    self.assertEqual(localElem, localElem2)
                    if elem.idx == localElem.idx:
                        self.assertTrue(localElem.isLocal())
                        self.assertFalse(elem == localElem)
                        self.assertFalse(elems[i:i+1] == elems[i:i+1].toLocal())

    def testDistCompPatchMembs(self):
        mesh = self.mesh3

        with mesh:
            TETGROUP1 = Compartment.Create(conductivity=1.23)
            TETGROUP2 = Compartment.Create()
            TETGROUP3 = Compartment.Create()
            TRIGROUP3 = Patch.Create(TETGROUP3)

            remSurfTris = TETGROUP1.surface & mesh.surface

            patch1 = Patch.Create(remSurfTris[:len(remSurfTris)//2], TETGROUP1)
            patch2 = Patch.Create(remSurfTris[len(remSurfTris)//2:], TETGROUP1)

            with self.assertRaises(ValueError):
                memb1 = Membrane.Create([patch1, patch2])

            memb2 = Membrane.Create([patch1], capacitance=4.56)
            memb3 = Membrane.Create([patch2], capacitance=7.89)

        self.assertEqual(set(TETGROUP1.tets), set(mesh.tetGroups['TETGROUP1']))
        self.assertEqual(set(TETGROUP2.tets), set(mesh.tetGroups['TETGROUP2']))
        self.assertEqual(set(TETGROUP3.tets), set(mesh.tetGroups['TETGROUP3']))
        self.assertEqual(set(TRIGROUP3.tris), set(mesh.triGroups['TRIGROUP3']))

        self.assertAlmostEqual(TETGROUP1.Vol + TETGROUP2.Vol + TETGROUP3.Vol, mesh.Vol)

        self.assertCountEqual(TETGROUP1.surface ^ TETGROUP2.surface ^ TETGROUP3.surface, mesh.surface)
        self.assertGreater(len(TETGROUP1.surface - mesh.surface), 0)
        self.assertGreater(len(TETGROUP2.surface - mesh.surface), 0)
        self.assertGreater(len(TETGROUP3.surface - mesh.surface), 0)

        with mesh.asLocal():
            tg1 = TETGROUP1.tets
            tg2 = TETGROUP2.tets
            tg3 = TETGROUP3.tets
            s1 = TETGROUP1.surface
            s2 = TETGROUP2.surface
            s3 = TETGROUP3.surface
            bb1 = TETGROUP1.bbox
            bb2 = TETGROUP2.bbox
            bb3 = TETGROUP3.bbox
            v1 = TETGROUP1.Vol
            v2 = TETGROUP2.Vol
            v3 = TETGROUP3.Vol

            trg3 = TRIGROUP3.tris
            p1 = patch1.tris
            p2 = patch2.tris
            tbb1 = TRIGROUP3.bbox
            tbb2 = patch1.bbox
            tbb3 = patch2.bbox
            a1 = TRIGROUP3.Area
            a2 = patch1.Area
            a3 = patch2.Area
        self.assertCountEqual(tg1, mesh.tetGroups['TETGROUP1'].toLocal())
        self.assertCountEqual(tg2, mesh.tetGroups['TETGROUP2'].toLocal())
        self.assertCountEqual(tg3, mesh.tetGroups['TETGROUP3'].toLocal())
        self.assertCountEqual(s1.combineWithOperator(operator.or_), TETGROUP1.surface)
        self.assertCountEqual(s2.combineWithOperator(operator.or_), TETGROUP2.surface)
        self.assertCountEqual(s3.combineWithOperator(operator.or_), TETGROUP3.surface)
        for locbb, expglobbb in [(bb1, TETGROUP1.bbox), (bb2, TETGROUP2.bbox), (bb3, TETGROUP3.bbox)]:
            globbbmin = self.allGather(locbb.min, lambda a,b: [min(x,y) for x,y in zip(a,b)])
            globbbmax = self.allGather(locbb.max, lambda a,b: [max(x,y) for x,y in zip(a,b)])
            self.assertEqual(globbbmin, list(expglobbb.min))
            self.assertEqual(globbbmax, list(expglobbb.max))
        self.assertEqual(self.allGather(v1), TETGROUP1.Vol)
        self.assertEqual(self.allGather(v2), TETGROUP2.Vol)
        self.assertEqual(self.allGather(v3), TETGROUP3.Vol)

        self.assertCountEqual(trg3, mesh.triGroups['TRIGROUP3'].toLocal())
        self.assertCountEqual(p1, patch1.tris.toLocal())
        self.assertCountEqual(p2, patch2.tris.toLocal())
        for locbb, expglobbb in [(tbb1, TRIGROUP3.bbox), (tbb2, patch1.bbox), (tbb3, patch2.bbox)]:
            globbbmin = self.allGather(locbb.min, lambda a,b: [min(x,y) for x,y in zip(a,b)])
            globbbmax = self.allGather(locbb.max, lambda a,b: [max(x,y) for x,y in zip(a,b)])
            self.assertEqual(globbbmin, list(expglobbb.min))
            self.assertEqual(globbbmax, list(expglobbb.max))
        self.assertEqual(self.allGather(a1), TRIGROUP3.Area)
        self.assertEqual(self.allGather(a2), patch1.Area)
        self.assertEqual(self.allGather(a3), patch2.Area)

        with mesh.asLocal(owned=False):
            with self.assertWarns(UserWarning):
                self.assertEqual(list(TETGROUP1.bbox.min), list(bb1.min))
            with self.assertWarns(UserWarning):
                self.assertEqual(list(TETGROUP1.bbox.max), list(bb1.max))

            with self.assertWarns(UserWarning):
                self.assertEqual(TETGROUP1.Vol, v1)
            with self.assertWarns(UserWarning):
                self.assertEqual(TETGROUP2.Vol, v2)
            with self.assertWarns(UserWarning):
                self.assertEqual(TETGROUP3.Vol, v3)

            with self.assertWarns(UserWarning):
                self.assertCountEqual(TETGROUP1.surface, s1)
            with self.assertWarns(UserWarning):
                self.assertCountEqual(TETGROUP2.surface, s2)
            with self.assertWarns(UserWarning):
                self.assertCountEqual(TETGROUP3.surface, s3)

            with self.assertWarns(UserWarning):
                self.assertEqual(TETGROUP1.Vol, v1)
            with self.assertWarns(UserWarning):
                self.assertEqual(TETGROUP2.Vol, v2)
            with self.assertWarns(UserWarning):
                self.assertEqual(TETGROUP3.Vol, v3)

            with self.assertWarns(UserWarning):
                self.assertEqual(list(TRIGROUP3.bbox.min), list(tbb1.min))
            with self.assertWarns(UserWarning):
                self.assertEqual(list(TRIGROUP3.bbox.max), list(tbb1.max))

            with self.assertWarns(UserWarning):
                self.assertEqual(TRIGROUP3.Area, a1)
            with self.assertWarns(UserWarning):
                self.assertEqual(patch1.Area, a2)
            with self.assertWarns(UserWarning):
                self.assertEqual(patch2.Area, a3)

        self.assertEqual(TETGROUP1.Conductivity, 1.23)

        self.assertEqual(set(memb2.tris), set(remSurfTris[:len(remSurfTris)//2]))
        self.assertEqual(set(memb3.tris), set(remSurfTris[len(remSurfTris)//2:]))
        self.assertEqual(memb2.Capacitance, 4.56)
        self.assertEqual(memb3.Capacitance, 7.89)

    def testDistCompFromLists(self):
        mesh = self.mesh3

        # Create a compartment from global list
        with mesh:
            TG1 = Compartment.Create(mesh.tetGroups['TETGROUP1'])
        self.assertCountEqual(TG1.tets, mesh.tetGroups['TETGROUP1'])
            
        # Check that compartments require local lists with ghost elements
        with self.assertRaises(Exception):
            with mesh.asLocal():
                tets3owned = mesh.tetGroups['TETGROUP3']
                TG2_ = Compartment.Create(mesh.tetGroups['TETGROUP2'])

        with mesh:
            with self.assertRaises(Exception):
                TG3_ = Compartment.Create(tets3owned)

        # Create a compartment from local list with ghost elements
        with mesh.asLocal(owned=False) as mesh:
            tets2 = mesh.tetGroups['TETGROUP2']
            TG2_1 = Compartment.Create(tets2[:len(tets2)//2])
            TG2_2 = Compartment.Create(tets2[len(tets2)//2:] & mesh.tets)
            tets3 = mesh.tetGroups['TETGROUP3']
        self.assertCountEqual(TG2_1.tets | TG2_2.tets, mesh.tetGroups['TETGROUP2'])

        with mesh:
            TG3 = Compartment.Create(tets3)
        self.assertCountEqual(TG3.tets, mesh.tetGroups['TETGROUP3'])

    def testDistPatchFromLists(self):
        mesh = self.mesh3

        with mesh:
            TETGROUP1 = Compartment.Create()
            TETGROUP2 = Compartment.Create()
            TETGROUP3 = Compartment.Create()

        # Create a patch from global lists
        with mesh:
            remSurfTris = TETGROUP1.surface & mesh.surface
            patch1 = Patch.Create(remSurfTris, TETGROUP1)
        self.assertCountEqual(patch1.tris, remSurfTris)

        # Check that patches require local lists with ghost bounds
        with self.assertRaises(Exception):
            with mesh.asLocal():
                tris2owned = TETGROUP2.surface & mesh.surface
                patch2_0 = Patch.Create(tris2owned, TETGROUP2)

        with mesh:
            with self.assertRaises(Exception):
                patch2_1 = Compartment.Create(tris2owned, TETGROUP2)

        # Create a patch from local list with ghost bounds
        tris2 = (TETGROUP2.surface & mesh.surface).toLocal(owned=False)
        tris3 = (TETGROUP3.surface & mesh.surface).toLocal(owned=False)
        with mesh.asLocal(owned=False) as mesh:
            patch2 = Patch.Create(tris2, TETGROUP2)
            self.assertCountEqual(patch2.tris, tris2)
        self.assertCountEqual(patch2.tris, tris2.combineWithOperator(operator.or_))

        with mesh:
            patch3 = Patch.Create(tris3, TETGROUP3)
        self.assertCountEqual(patch3.tris, tris3.combineWithOperator(operator.or_))

        # Create inner patch
        trisInner = TETGROUP2.surface & TETGROUP3.surface
        with mesh.asLocal():
            with self.assertRaises(Exception):
                patchInner_1 = Patch.Create(trisInner.toLocal(), TETGROUP2, TETGROUP3)
            patchInner = Patch.Create(trisInner.toLocal(owned=False), TETGROUP2, TETGROUP3)
            self.assertCountEqual(patchInner.tris, trisInner.toLocal())
        with mesh.asLocal(owned=False):
            self.assertCountEqual(patchInner.tris, trisInner.toLocal(owned=False))
        self.assertCountEqual(patchInner.tris, trisInner)

    def testSplitMeshElemNumbering_n4(self):
        # Here we load all parts, there should not be any jump in mesh indexes
        splitMesh = DistMesh(os.path.join(FILEDIR, '../../../../mesh/CaBurst/branch_split'))

        # Check that the element ids are between 0 and N - 1
        self.assertCountEqual(splitMesh.tets.indices, splitMesh.stepsMesh.getAllTetIndices())
        self.assertCountEqual(splitMesh.tris.indices, splitMesh.stepsMesh.getAllTriIndices())
        self.assertCountEqual(splitMesh.verts.indices, splitMesh.stepsMesh.getAllVertIndices())

    def testSplitMeshElemNumberingJumps_n2(self):
        # We load on purpose only two parts out of 4 to be sure there will be jumps in element numbers
        splitMesh = DistMesh(os.path.join(FILEDIR, '../../../../mesh/CaBurst/branch_split'))

        # Check that the element ids are between 0 and N - 1
        self.assertCountEqual(splitMesh.tets.indices, splitMesh.stepsMesh.getAllTetIndices())
        self.assertCountEqual(splitMesh.tris.indices, splitMesh.stepsMesh.getAllTriIndices())
        self.assertCountEqual(splitMesh.verts.indices, splitMesh.stepsMesh.getAllVertIndices())

        # Check that they are also correctly ordered
        self.assertEqual(list(splitMesh.tets.indices), splitMesh.stepsMesh.getAllTetIndices())
        self.assertEqual(list(splitMesh.tris.indices), splitMesh.stepsMesh.getAllTriIndices())
        self.assertEqual(list(splitMesh.verts.indices), splitMesh.stepsMesh.getAllVertIndices())

    def testgetTriTetNeighbs(self):
        # we avoid the vertex that is contexted by multiple tets as the answer is not unique
        centerTet = self.mesh.tets[1e-10, 1e-10, 1e-10]
        centerTri = centerTet.faces[0]
        tetNeighbs = centerTri.tetNeighbs
        self.assertEqual(len(tetNeighbs), 2)
        self.assertTrue(centerTet in tetNeighbs)

        surfTet = self.mesh.surface[0]
        self.assertEqual(len(surfTet.tetNeighbs), 1)


del test_tetMesh.tetMeshTests.testVTKLoading
del test_tetMesh.tetMeshTests.testGmshLoading
del test_tetMesh.tetMeshTests.testTetGenLoading
del test_tetMesh.tetMeshTests.testAbaqusLoading
del test_tetMesh.tetMeshTests.testLoading


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(distTetMeshTests, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
