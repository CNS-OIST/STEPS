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

""" Unit tests for tetmesh class and related methods."""

import unittest
import tempfile
import os

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.utils import *

from steps.API_1.geom import UNKNOWN_TET, UNKNOWN_TRI

FILEDIR = os.path.dirname(os.path.abspath(__file__))

class tetMeshTests(unittest.TestCase):
    """Test tetmesh loading and geometrical elements handling."""
    def setUp(self):
        self.setUpMeshes()

        self.model = Model()
        with self.model:
            vsys = VolumeSystem.Create()
            ssys = SurfaceSystem.Create()

            SA, SB, SC = Species.Create()
            r = ReactionManager()
            with vsys:
                SA + SB <r[0]> SC
                r[0].K = 1,1
            with ssys:
                SA.s + SB.s <r[0]> SC.s
                r[0].K = 1,1

            self.vsys, self.ssys, self.SA, self.SB, self.SC = vsys, ssys, SA, SB, SC

    def setUpMeshes(self):
        self._distMesh = False
        self.refMesh = TetMesh.LoadAbaqus(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets.inp'), 1e-6)

        self.mesh = TetMesh.Load(os.path.join(FILEDIR, 'meshes', 'cyl_len10_diam1'))
        self.mesh2 = TetMesh.Load(os.path.join(FILEDIR, 'meshes', 'cyl_len10_diam1'))

    def getSimulation(self):
        rng = RNG('mt19937', 512, 7233)
        return Simulation('Tetexact', self.model, self.mesh, rng)

    def checkIdenticalMeshes(self, mesh1, mesh2):
        self.assertEqual(len(mesh1.tets), len(mesh2.tets))
        self.assertEqual(len(mesh1.tris), len(mesh2.tris))

        self.assertListEqual([tet.Vol for tet in mesh1.tets], [tet.Vol for tet in mesh2.tets])
        self.assertListEqual([tri.Area for tri in mesh1.tris], [tri.Area for tri in mesh2.tris])

    def testVTKLoading(self):
        mesh = TetMesh.LoadVTK(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets.vtk'), 1e-6)
        self.checkIdenticalMeshes(self.refMesh, mesh)

    def testGmshLoading(self):
        mesh = TetMesh.LoadGmsh(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets.msh'), 1e-6)
        self.checkIdenticalMeshes(self.refMesh, mesh)

    def testTetGenLoading(self):
        mesh = TetMesh.LoadTetGen(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets'), 1e-6)
        self.checkIdenticalMeshes(self.refMesh, mesh)

    def testAbaqusLoading(self):
        with self.assertRaises(TypeError):
            mesh = TetMesh.LoadAbaqus(42, 1e-6)

        # Abaqus2
        tetFile = os.path.join(FILEDIR, 'meshes', 'mesh_tet.inp')
        triFile = os.path.join(FILEDIR, 'meshes', 'mesh_tri.inp')
        mesh = TetMesh.LoadAbaqus((tetFile, triFile), 1e-6)

        self.assertEqual(len(mesh.tets), 535)

    def testGroupBlocksLoading(self):
        mesh = TetMesh.LoadAbaqus(os.path.join(FILEDIR, 'meshes', 'sphere_mesh.inp'), 1e-6)

        self.assertSetEqual(set(mesh.vertGroups.keys()), set(['NODEGROUP1', 'NODEGROUP2']))
        self.assertSetEqual(set(mesh.triGroups.keys()), set(['TRIGROUP1', 'TRIGROUP2']))
        self.assertSetEqual(set(mesh.tetGroups.keys()), set(['TETGROUP1', 'TETGROUP2']))

        tets = set(mesh.tetGroups['TETGROUP1'] + mesh.tetGroups['TETGROUP2'])
        self.assertSetEqual(tets, set(mesh.tets))

        tris = set(mesh.triGroups['TRIGROUP1'] + mesh.triGroups['TRIGROUP2'])
        trisref = set(mesh.tetGroups['TETGROUP1'].surface)
        self.assertSetEqual(tris, trisref)

        nodes = set(mesh.vertGroups['NODEGROUP1'] + mesh.vertGroups['NODEGROUP2'])
        nodesref = set(mesh.tets.verts)
        self.assertSetEqual(nodes, nodesref)

    def testLoading(self):
        """Test loading of tetmesh from other file formats."""
        # Prevent creation of bare tetmesh
        with self.assertRaises(Exception):
            TetMesh()

        # Abaqus
        mesh = TetMesh.LoadAbaqus(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets.inp'), 1e-6)
        self.assertEqual(len(mesh.tets), 1400)

        # Loading saved compartments and patches
        with mesh:
            n = len(mesh.tets)
            c1tets, c2tets = mesh.tets[:n//2], mesh.tets[n//2:]
            comp1 = Compartment.Create(c1tets, self.model.ssys)
            comp2 = Compartment.Create(c2tets, self.model.ssys)
            p1tris = c1tets.surface & c2tets.surface
            patch1 = Patch.Create(p1tris, comp1, comp2, self.model.ssys)
            roi1 = ROI.Create(mesh.tets[n//4:-n//4])

        with tempfile.TemporaryDirectory() as tmpdir:
            _, meshPath = tempfile.mkstemp(dir=tmpdir, prefix='meshName')
            mesh.Save(meshPath)

            #Loading the mesh and the compartments / patches
            newMesh2 = TetMesh.Load(meshPath)
            self.assertEqual(len(newMesh2.tets), len(mesh.tets))
            self.assertEqual(len(newMesh2.tris), len(mesh.tris))
            self.assertEqual(set(tet.idx for tet in newMesh2.comp1.tets), set(tet.idx for tet in mesh.comp1.tets))
            self.assertEqual(set(tet.idx for tet in newMesh2.comp2.tets), set(tet.idx for tet in mesh.comp2.tets))
            self.assertEqual(set(tet.idx for tet in newMesh2.patch1.tris), set(tet.idx for tet in mesh.patch1.tris))
            self.assertEqual(set(vsys.name for vsys in newMesh2.comp1.systems), set(vsys.name for vsys in mesh.comp1.systems))
            self.assertEqual(set(vsys.name for vsys in newMesh2.comp2.systems), set(vsys.name for vsys in mesh.comp2.systems))
            self.assertEqual(set(vsys.name for vsys in newMesh2.patch1.systems), set(vsys.name for vsys in mesh.patch1.systems))
            self.assertEqual(set(tet.idx for tet in newMesh2.roi1.tets), set(tet.idx for tet in mesh.roi1.tets))


    def testBasicAccess(self):
        """Test access to geometrical elements directly from the mesh."""
        stepsMesh = self.mesh.stepsMesh
        # Iteration
        self.assertEqual(tuple(tet.idx for tet in self.mesh.tets), tuple(range(stepsMesh.countTets())))
        self.assertEqual(tuple(tet.idx for tet in self.mesh.tris), tuple(range(stepsMesh.countTris())))
        if not self._distMesh:
            self.assertEqual(tuple(tet.idx for tet in self.mesh.bars), tuple(range(stepsMesh.countBars())))
        self.assertEqual(tuple(tet.idx for tet in self.mesh.verts), tuple(range(stepsMesh.countVertices())))

        # Access with []
        n = stepsMesh.countTets()
        self.assertEqual(self.mesh.tets[n // 2].idx, n // 2)
        with self.assertRaises(IndexError):
            self.mesh.tets[2*n]
        self.assertEqual(self.mesh.tets[0, 0, 0].idx, stepsMesh.findTetByPoint([0, 0, 0]))
        with self.assertRaises(KeyError):
            self.mesh.tets[10, 10, 10]

        n = stepsMesh.countTris()
        self.assertEqual(self.mesh.tris[n // 2].idx, n // 2)
        with self.assertRaises(IndexError):
            self.mesh.tris[2*n]

        if not self._distMesh:
            n = stepsMesh.countBars()
            self.assertEqual(self.mesh.bars[n // 2].idx, n // 2)
            with self.assertRaises(IndexError):
                self.mesh.bars[2*n]

        n = stepsMesh.countVertices()
        self.assertEqual(self.mesh.verts[n // 2].idx, n // 2)
        with self.assertRaises(IndexError):
            self.mesh.verts[2*n]

        # Mesh properties
        self.mesh.Vol

    def testElemMethods(self):
        """Test geometrical elements method to access faces from tets, bars from triangles, etc."""
        stepsMesh = self.mesh.stepsMesh

        # Tet
        tet = self.mesh.tets[0, 0, 0]
        tet2 = TetReference([0, 0, 0], mesh=self.mesh)
        self.assertEqual(tet, tet2)
        self.assertEqual(tuple(tet.center), tuple(stepsMesh.getTetBarycenter(tet.idx)))
        self.assertEqual(tuple(tri.idx for tri in tet.faces), tuple(stepsMesh.getTetTriNeighb(tet.idx)))
        self.assertEqual(tet.Vol, stepsMesh.getTetVol(tet.idx))
        if not self._distMesh:
            self.assertEqual(tet.qualityRER, stepsMesh.getTetQualityRER(tet.idx))
        self.assertEqual(tuple(tet2.idx for tet2 in tet.neighbs), tuple(ind for ind in stepsMesh.getTetTetNeighb(tet.idx) if ind != UNKNOWN_TET))
        self.assertEqual(tuple(v.idx for v in tet.verts), tuple(stepsMesh.getTet(tet.idx)))
        self.assertEqual(tet.toList(), TetList([tet], mesh=self.mesh))
        with self.assertRaises(IndexError):
            tet.faces[4]
        with self.assertRaises(IndexError):
            tet.neighbs[4]
        with self.assertRaises(IndexError):
            tet.verts[4]
        self.assertTrue(tet.containsPoint([0, 0, 0]))
        self.assertFalse(tet.containsPoint([1, 0, 0]))

        # Tri
        tri = tet.faces[0]
        self.assertEqual(tuple(tri.center), tuple(stepsMesh.getTriBarycenter(tri.idx)))
        self.assertEqual(tri.Area, stepsMesh.getTriArea(tri.idx))
        self.assertEqual(tuple(tet2.idx for tet2 in tri.tetNeighbs), tuple(ind for ind in stepsMesh.getTriTetNeighb(tri.idx) if ind != UNKNOWN_TET))
        self.assertEqual(tuple(v.idx for v in tri.verts), tuple(stepsMesh.getTri(tri.idx)))
        self.assertEqual(tri.toList(), TriList([tri], mesh=self.mesh))
        tetNeighbTris = TriList([], mesh=self.mesh)
        neighbs = tri.tetNeighbs
        for t in neighbs:
            tetNeighbTris |= t.faces
        tetNeighbTris -= tri.toList()
        self.assertGreaterEqual(len(tri.triNeighbs), len(tetNeighbTris))
        self.assertEqual(len(tetNeighbTris - tri.triNeighbs), 0)
        if not self._distMesh:
            self.assertEqual(tuple(b.idx for b in tri.bars), tuple(stepsMesh.getTriBars(tri.idx)))
            self.assertEqual(tuple(tri.norm), tuple(stepsMesh.getTriNorm(tri.idx)))
            with self.assertRaises(IndexError):
                tri.bars[3]
        with self.assertRaises(IndexError):
            tri.tetNeighbs[2]
        with self.assertRaises(IndexError):
            tri.verts[3]

        # Bar
        if not self._distMesh:
            bar = tri.bars[0]
            v1, v2 = bar.verts
            self.assertEqual(tuple(bar.center), tuple((v1 + v2)/2))
            self.assertEqual(tuple(v.idx for v in bar.verts), tuple(stepsMesh.getBar(bar.idx)))
            self.assertEqual(bar.toList(), BarList([bar], mesh=self.mesh))
            with self.assertRaises(IndexError):
                bar.verts[2]

            # Vert
            v = bar.verts[0]
            self.assertEqual((v.x, v.y, v.z), tuple(stepsMesh.getVertex(v.idx)))
            self.assertEqual((v[0], v[1], v[2]), tuple(stepsMesh.getVertex(v.idx)))
            self.assertEqual(v.toList(), VertList([v], mesh=self.mesh))
            with self.assertRaises(IndexError):
                v[3]

        # Build from other references
        tet1 = self.mesh.tets[0]
        tet2 = TetReference(tet1)
        tet3 = TetReference(0, mesh=self.mesh2)
        self.assertEqual(tet1, tet2)

        with self.assertRaises(Exception):
            TetReference(tet1.idx)

        with self.assertRaises(Exception):
            TetReference(tet3, mesh=self.mesh)

    def testElemLists(self):
        """
        Test the different ways to create, access, and combine lists of geometrical elements.
        """
        def elemTest(lstCls, allElems, nbAppend = 10):
            otherCls = TriList if lstCls == TetList else TetList
            # Test reference creation
            self.assertEqual(lstCls._refCls(0, self.mesh), allElems[0])
            with self.assertRaises(TypeError):
                lstCls._refCls(10.0, self.mesh)
            with self.assertRaises(TypeError):
                lstCls._refCls(0, 10)

            # test creating from scratch
            with self.assertRaises(TypeError):
                lstCls(None, 'mesh')
            lst = lstCls(None, self.mesh)
            self.assertEqual(len(lst), 0)
            with self.assertRaises(IndexError):
                lst[0]
            lst2 = lstCls([2 * i + 1 for i in range(nbAppend)], self.mesh)
            with self.assertRaises(TypeError):
                lst2 = lstCls(10, self.mesh)
            with self.assertRaises(TypeError):
                lst2 = lstCls([0, 1, 5.2], self.mesh)
            lst3 = lstCls(range(nbAppend), self.mesh)
            lst4 = lstCls(list(range(nbAppend)), self.mesh)
            lst5 = lstCls(lst4, self.mesh)
            lst6 = lstCls([e for e in allElems[0:nbAppend]], self.mesh)
            lst7 = lstCls(range(nbAppend), self.mesh)
            lst8 = lstCls(range(nbAppend), self.mesh2)
            lst9 = lstCls(mesh=self.mesh)

            # test creation from generator
            lg1 = lstCls((i for i in range(nbAppend)), self.mesh)
            self.assertEqual(len(lg1), nbAppend)
            self.assertEqual(lg1.indices, list(range(nbAppend)))

            # test appending
            with self.assertRaises(TypeError):
                lst.append(0)
            lst.append(allElems[0])
            self.assertEqual(len(lst), 1)
            self.assertEqual(lst[0].idx, 0)
            for i in range(nbAppend):
                lst.append(allElems[2 * i + 1])
            lst7.append(allElems[nbAppend])

            # test removing
            with self.assertRaises(TypeError):
                lst9.remove(0)
            with self.assertRaises(ValueError):
                lst9.remove(allElems[-1])
            lst9.append(allElems[-1])
            self.assertEqual(len(lst9), 1)
            lst9.remove(allElems[-1])
            self.assertEqual(len(lst9), 0)

            # test arbitrary position accessing
            self.assertEqual(lst[nbAppend].idx, 2 * nbAppend - 1)
            self.assertEqual(lst[1].idx, 1)

            # Test accesing with -1
            self.assertEqual(lst[-1].idx, 2 * nbAppend - 1)
            self.assertEqual(lst[-2].idx, 2 * nbAppend - 3)
            with self.assertRaises(IndexError):
                lst[-10*nbAppend]

            # Test illegal accessing
            with self.assertRaises(TypeError):
                lst[None]

            # test equality of references
            self.assertEqual(lst[0], allElems[0])
            self.assertEqual(lst[-1], allElems[2 * nbAppend - 1])
            self.assertNotEqual(lst[0], lst[1])
            self.assertNotEqual(lst[1], lst[-1])

            # test presence of reference in lists
            self.assertIn(allElems[0], lst)
            self.assertIn(allElems[2 * nbAppend - 1], lst)
            self.assertIn(allElems[1], lst)
            self.assertNotIn(0, lst)
            self.assertNotIn(otherCls._refCls(0, self.mesh), lst3)

            # test equality of lists
            self.assertEqual(lst[1:], lst2)
            self.assertEqual(lst3, lst4)
            self.assertEqual(lst4, lst5)
            self.assertEqual(lst5, lst6)
            lst4.append(allElems[nbAppend])
            self.assertNotEqual(lst3, lst4)
            self.assertNotEqual(lst3, lst2)
            self.assertNotEqual(lst8, lst3)

            # test combinations
            lA = lstCls(range(0, 6), self.mesh)
            lB = lstCls(range(3, 9), self.mesh)
            lAiB = lstCls(range(3, 6), self.mesh)
            lAuB = lstCls(range(9), self.mesh)
            lAmB = lstCls(range(0, 3), self.mesh)
            lApB = lstCls(list(range(0, 6)) + list(range(3, 9)), self.mesh)
            lAxB = lstCls([0, 1, 2, 6, 7, 8], self.mesh)

            self.assertEqual(lA & lB, lAiB)
            self.assertEqual(lA | lB, lAuB)
            self.assertEqual(lA - lB, lAmB)
            self.assertEqual(lA + lB, lApB)
            self.assertEqual(lA ^ lB, lAxB)

            lC = otherCls(range(3, 9), self.mesh)
            with self.assertRaises(TypeError):
                lA & lC
            with self.assertRaises(TypeError):
                lA | lC
            with self.assertRaises(TypeError):
                lA - lC
            with self.assertRaises(TypeError):
                lA + lC
            with self.assertRaises(TypeError):
                lA ^ lC
            with self.assertRaises(Exception):
                lst8 + lA

            # check immutability of mesh lists
            with self.assertRaises(Exception):
                allElems.append(allElems[0])

        nbAppend = 10
        elemTest(TetList, self.mesh.tets, nbAppend)
        elemTest(TriList, self.mesh.tris, nbAppend)
        if not self._distMesh:
            elemTest(BarList, self.mesh.bars, nbAppend)
        elemTest(VertList, self.mesh.verts, nbAppend)

        # check non-assignability of mesh lists
        with self.assertRaises(Exception):
            self.mesh.tets += self.mesh.tets
        with self.assertRaises(Exception):
            self.mesh.tris += self.mesh.tris
        with self.assertRaises(Exception):
            self.mesh.bars += self.mesh.bars
        with self.assertRaises(Exception):
            self.mesh.verts += self.mesh.verts

        # Tests specific to TetList
        lst1 = TetList(range(nbAppend), self.mesh)
        lst2 = TetList(None, self.mesh)
        lst2.append(self.mesh.tets[0, 0, 0])
        lst1 -= lst2

        with self.assertRaises(KeyError):
            lst1[0, 0, 0]
        self.assertNotIn((0, 0, 0), lst1)
        self.assertIn((0, 0, 0), lst2)
        self.assertNotIn(Point(0, 0, 0), lst1)
        self.assertIn(Point(0, 0, 0), lst2)

        tetSurf = TriList(sorted(tri.idx for tri in self.mesh.tets.surface), self.mesh)
        meshSurf = TriList(sorted(tri.idx for tri in self.mesh.surface), self.mesh)
        self.assertTrue(len(tetSurf) > 0)
        self.assertEqual(tetSurf, meshSurf)

        ## .tris
        lst1t = lst1.tris
        checked = [False] * len(lst1t)
        for tet in lst1:
            for tri in tet.faces:
                self.assertTrue(tri in lst1t)
                checked[lst1t.index(tri)] = True
        self.assertTrue(all(checked))

        ## .bars
        if not self._distMesh:
            lst1b = lst1.bars
            checked = [False] * len(lst1b)
            for tet in lst1:
                for tri in tet.faces:
                    for bar in tri.bars:
                        self.assertTrue(bar in lst1b)
                        checked[lst1b.index(bar)] = True
            self.assertTrue(all(checked))

        ## .verts
        lst1v = lst1.verts
        checked = [False] * len(lst1v)
        for tet in lst1:
            for tri in tet.faces:
                if not self._distMesh:
                    for bar in tri.bars:
                        for vert in bar.verts:
                            self.assertTrue(vert in lst1v)
                            checked[lst1v.index(vert)] = True
                else:
                    for vert in tri.verts:
                        self.assertTrue(vert in lst1v)
                        checked[lst1v.index(vert)] = True
        self.assertTrue(all(checked))

        vol1 = sum(tet.Vol for tet in lst1)
        self.assertAlmostEqual(lst1.Vol, vol1)
        self.assertAlmostEqual((lst1 + lst1).Vol, vol1)

        ## .dilate() and .erode()
        lst1 = TetList([self.mesh.tets[0, 0, 0]])
        lst2 = TetList(lst1)
        lst2.dilate(0)
        self.assertEqual(lst1, lst2)
        lst2.dilate(1)
        self.assertNotEqual(lst1, lst2)
        lst2.erode(1)
        if len(lst2) > 0:
            self.assertEqual(lst1, lst2)

        lst2.erode(len(lst2))
        self.assertEqual(len(lst2), 0)

        lst2 = TetList(lst1)
        lst2.dilate(len(self.mesh.tets))
        self.assertEqual(set(lst2), set(self.mesh.tets))

        lst2 = TetList(self.mesh.tets)
        lst2.erode(1)
        lst2 = self.mesh.tets - lst2

        self.assertEqual(set(lst2.surface & self.mesh.surface), set(self.mesh.surface))

        with self.assertRaises(Exception):
            self.mesh.tets.dilate(1)

        with self.assertRaises(Exception):
            self.mesh.tets.erode(1)

        # Tests specific to TriList
        if not self._distMesh:
            self.assertEqual(len(tetSurf.edges), 0)
            self.assertTrue(len(tetSurf[0:nbAppend].edges) > 0)

        lst1 = TriList(range(nbAppend), self.mesh)

        ## .bars
        if not self._distMesh:
            lst1b = lst1.bars
            checked = [False] * len(lst1b)
            for tri in lst1:
                for bar in tri.bars:
                    self.assertTrue(bar in lst1b)
                    checked[lst1b.index(bar)] = True
            self.assertTrue(all(checked))

        ## .verts
        lst1v = lst1.verts
        checked = [False] * len(lst1v)
        for tri in lst1:
            if not self._distMesh:
                for bar in tri.bars:
                    for vert in bar.verts:
                        self.assertTrue(vert in lst1v)
                        checked[lst1v.index(vert)] = True
            else:
                for vert in tri.verts:
                    self.assertTrue(vert in lst1v)
                    checked[lst1v.index(vert)] = True
        self.assertTrue(all(checked))

        area1 = sum(tri.Area for tri in lst1)
        self.assertAlmostEqual(lst1.Area, area1)
        self.assertAlmostEqual((lst1 + lst1).Area, area1)

        # Tests specific to BarList

        if not self._distMesh:
            lst1 = BarList(range(nbAppend), self.mesh)

            ## .verts
            lst1v = lst1.verts
            checked = [False] * len(lst1v)
            for bar in lst1:
                for vert in bar.verts:
                    self.assertTrue(vert in lst1v)
                    checked[lst1v.index(vert)] = True
            self.assertTrue(all(checked))

    def testCompartments(self):
        """Test the creation of compartments and data acces to/from tetrahedrons."""
        # compartment creation
        nbComps = 6
        notAttrInd = nbComps - 1
        nbTetPerComp = len(self.mesh.tets) // nbComps
        tetLsts = [inds for inds in zip(*([iter(self.mesh.tets)]*nbTetPerComp))]
        
        # Throw exception if declaring outside a mesh
        with self.assertRaises(Exception):
            comp0 = Compartment.Create(tetLsts[0])
        # Normal creation
        with self.mesh:
            with self.assertRaises(Exception):
                compA = Compartment.Create(self.vsys)

            comp1 = Compartment.Create(tetLsts[1])
            comp1.addSystem(self.vsys)

            comp2 = Compartment.Create(tetLsts[2], self.vsys)

            comp3 = Compartment.Create(TetList(tetLsts[3], self.mesh))

            with self.assertRaises(Exception):
                Compartment(TriList(range(5), self.mesh))
            if not self._distMesh:
                with self.assertRaises(Exception):
                    comp4 = Compartment.Create(tetLsts[4], self.vsys, 1)
                with self.assertRaises(Exception):
                    Compartment()
            with self.assertRaises(Exception):
                Compartment(5)
            with self.assertRaises(TypeError):
                comp1.addSystem(1)

            if self._distMesh:
                # DistTetOpSplit requires that all tetrahedrons are associated with a compartment
                notAttrComp = Compartment.Create(tetLsts[0] + tetLsts[4] + tetLsts[notAttrInd])

        # checking that tetrahedrons have their comp thing working
        self.assertEqual(tetLsts[1][0].comp, comp1)
        self.assertNotEqual(tetLsts[1][0].comp, comp2)
        if not self._distMesh:
            self.assertIsNone(tetLsts[notAttrInd][0].comp)

        # accessing tets from a comp
        self.assertCountEqual(comp1.tets, TetList(tetLsts[1], self.mesh))
        self.assertCountEqual(comp2.tets, TetList(tetLsts[2], self.mesh))
        self.assertEqual(len(set(comp2.tets.indices) & set(tetLsts[3])), 0)

        # Declare simulation to allow species access from patch
        sim = self.getSimulation()

        # accessing species from a comp
        self.assertEqual(comp1.SA, self.SA)
        self.assertEqual(comp1.SB, self.SB)
        with self.assertRaises(Exception):
            comp1.SD

        # accessing species from a tet
        self.assertEqual(tetLsts[1][0].SA, self.SA)
        self.assertEqual(tetLsts[1][0].SB, self.SB)
        with self.assertRaises(Exception):
            tetLsts[1][0].SD
        if not self._distMesh:
            with self.assertRaises(Exception):
                tetLsts[notAttrInd][0].SA

        # bbox
        c1b = comp1.bbox
        mb = self.mesh.bbox
        self.assertTrue(mb.min.x <= c1b.min.x)
        self.assertTrue(mb.min.y <= c1b.min.y)
        self.assertTrue(mb.min.z <= c1b.min.z)
        self.assertTrue(mb.max.x >= c1b.max.x)
        self.assertTrue(mb.max.y >= c1b.max.y)
        self.assertTrue(mb.max.z >= c1b.max.z)

        comp1.Vol

    def testPatches(self):
        """Test the creation of patches and data acces to/from triangles."""
        # prepare compartments
        # c1lst = TetList([tet for tet in self.mesh.tets if tet.center.x >= self.mesh.bbox.center.x], self.mesh)
        # TODO Temporary fix to avoid the "RuntimeError: SReac : the outer compartment element 398 is not
        #      in the same process as the patch patch3" errors.
        c1lst = TetList([tet for tet in self.mesh.tets if tet.center.x >= self.mesh.bbox.center.x * 0.3], self.mesh)
        c2lst = self.mesh.tets - c1lst
        trilst = c1lst.surface & c2lst.surface

        with self.mesh:
            comp1 = Compartment.Create(c1lst, self.vsys)
            comp2 = Compartment.Create(c2lst, self.vsys)

        # compartment creation
        nbComps = 6
        notAttrInd = nbComps - 1
        nbTriPerComp = len(trilst) // nbComps
        triLsts = [inds for inds in zip(*([iter(trilst)]*nbTriPerComp))]
        
        # Throw exception if declaring outside a mesh
        with self.assertRaises(Exception):
            patch0 = Patch.Create(triLsts[0], comp1, comp2, self.ssys)
        # Normal creation
        with self.mesh:
            with self.assertRaises(Exception):
                patchA = Patch.Create(comp1, comp2)

            patch1 = Patch.Create(triLsts[1], comp1, comp2)
            patch1.addSystem(self.ssys)

            patch2 = Patch.Create(triLsts[2], comp1, comp2, self.ssys)

            patch3 = Patch.Create(TriList(triLsts[3], self.mesh), comp1, comp2)

            if not self._distMesh:
                # Only check whith tetexact since otherwise the exception is only raised in a single rank
                with self.assertRaises(Exception):
                    Patch(TriList(triLsts[2][0:5], self.mesh), comp1, comp2)
            # TODO Not sure yet if setting the area of a tetPatch should raise an exception, the C++ side doesnt raise a NotImplError
            # with self.assertRaises(Exception):
                # patch4 = Patch.Create(triLsts[4], comp1, comp2, self.ssys, 1)
            with self.assertRaises(Exception):
                Patch()
            with self.assertRaises(Exception):
                Patch(5)
            with self.assertRaises(TypeError):
                patch1.addSystem(1)

        # checking that tetrahedrons have their patch thing working
        self.assertEqual(triLsts[1][0].patch, patch1)
        self.assertNotEqual(triLsts[1][0].patch, patch2)
        self.assertIsNone(triLsts[notAttrInd][0].patch)

        # accessing tets from a patch
        self.assertEqual(patch1.tris, TriList(triLsts[1], self.mesh))
        self.assertEqual(patch2.tris, TriList(triLsts[2], self.mesh))
        self.assertNotEqual(patch2.tris, TriList(triLsts[3], self.mesh))

        # Triangle neighbors in patch
        for tri in patch1.tris:
            self.assertCountEqual(tri.patchTriNeighbs, tri.triNeighbs & patch1.tris)
        self.assertEqual(self.mesh.tris[0].patchTriNeighbs, TriList([], self.mesh))

        # edges of a patch
        if not self._distMesh:
            self.assertEqual(set(patch1.edges), set(TriList(triLsts[1], self.mesh).edges))

        # Declare simulation to allow species access from patch
        sim = self.getSimulation()

        # accessing species from a patch
        self.assertEqual(patch1.SA, self.SA)
        self.assertEqual(patch1.SB, self.SB)
        with self.assertRaises(Exception):
            patch1.SD

        # accessing species from a tet
        self.assertEqual(triLsts[1][0].SA, self.SA)
        self.assertEqual(triLsts[1][0].SB, self.SB)
        with self.assertRaises(Exception):
            triLsts[1][0].SD
        with self.assertRaises(Exception):
            triLsts[notAttrInd][0].SA

        # bbox
        c1b = patch1.bbox
        mb = self.mesh.bbox
        self.assertTrue(mb.min.x <= c1b.min.x)
        self.assertTrue(mb.min.y <= c1b.min.y)
        self.assertTrue(mb.min.z <= c1b.min.z)
        self.assertTrue(mb.max.x >= c1b.max.x)
        self.assertTrue(mb.max.y >= c1b.max.y)
        self.assertTrue(mb.max.z >= c1b.max.z)

        patch1.Area


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(tetMeshTests, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
