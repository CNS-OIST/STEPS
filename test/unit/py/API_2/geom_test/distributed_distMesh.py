####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

""" Unit tests for distmesh class and related methods."""

import functools
import unittest
import tempfile
import os

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

        self.elemListNames = ['tets', 'tris', 'verts']

    def getSimulation(self):
        rng = RNG('mt19937', 512, 7233)
        return Simulation('DistTetOpSplit', self.model, self.mesh, rng)

    def gather(self, lst, func=lambda a, b: a + b):
        lists = self.mesh3._comm.gather(lst, root=0)
        if lists is not None:
            return functools.reduce(func, lists)
        else:
            return None

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
                self.assertLessEqual(l, g)
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

        # Vol
        vol = mesh.Vol
        with mesh.asLocal():
            localVol = mesh.Vol
            self.assertLessEqual(localVol, vol)
            totVol = self.gather(localVol)
            if totVol is not None:
                self.assertEqual(totVol, vol)
        self.assertEqual(vol, mesh.Vol)

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

        self.assertEqual(TETGROUP1.Conductivity, 1.23)

        self.assertEqual(set(memb2.tris), set(remSurfTris[:len(remSurfTris)//2]))
        self.assertEqual(set(memb3.tris), set(remSurfTris[len(remSurfTris)//2:]))
        self.assertEqual(memb2.Capacitance, 4.56)
        self.assertEqual(memb3.Capacitance, 7.89)

    def testSplitMeshElemNumbering(self):
        splitMesh = DistMesh(os.path.join(FILEDIR, '../../../../mesh/3tets_2patches_2comp_split2/3tets_2patches_2comp'))

        # Check that the element ids work as expected
        self.assertEqual(set(splitMesh.tets.indices), set(splitMesh.stepsMesh.getAllTetIndices()))
        self.assertEqual(set(splitMesh.tris.indices), set(splitMesh.stepsMesh.getAllTriIndices()))
        self.assertEqual(set(splitMesh.verts.indices), set(splitMesh.stepsMesh.getAllVertIndices()))


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
