####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

""" Unit tests for DistTetOpSplit.setVertIClamp(...)"""

import os
import unittest

import numpy as np

import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class DistTetopsplitEFieldClamps(unittest.TestCase):
    """Test that Efield clamping (for current and potential) works as expected with DistTetopsplit"""

    def setUp(self):
        self.Iinj = 1e-15
        self.capacitance = 0.01
        self.setUpModel()
        self.setUpMeshes()
        self.setUpSimulation()

    def setUpModel(self):
        """Set up model with reactions and diffusions"""
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA = Species.Create()
            vsys = VolumeSystem.Create()
            with vsys:
                SA >r[1]> None
                r[1].K = 1
        self.mdl = mdl

    def setUpMeshes(self):
        """ Create mesh and assign to it the VolumeSystem"""
        self.mesh = DistMesh(os.path.join(MESH_DIR, 'box.msh'), scale=1e-6)
        with self.mesh:
            comp = Compartment.Create(self.mesh.tets, self.mdl.vsys)
            surf = Patch.Create(self.mesh.surface, comp, None)
            membrane = Membrane.Create([surf], capacitance=self.capacitance)
            comp.Conductivity = 1

        self.area = surf.Area
        self.injVerts = self.mesh.surface.verts
        self.injTris = self.mesh.surface
        self.nbVerts = len(self.injVerts)
        self.nbSurfTris = len(self.injTris)

    def setUpSimulation(self):
        """ Instantiate main simulator object"""
        rng = RNG('mt19937', 512, 7233)
        self.sim = Simulation('DistTetOpSplit', self.mdl, self.mesh, rng)
        self.sim.EfieldDT = 1e-3

    def test_disttetopsplit_setVertIClamp_n2(self):
        """ Test after ENDTIME that the potential corresponds to what is expected given the current that should be injected"""
        ENDTIME = 0.1

        self.sim.newRun()

        self.sim.membrane.Potential = 0
        self.sim.VERTS(self.injVerts).IClamp = self.Iinj

        self.sim.run(ENDTIME)

        avgPot = sum(self.sim.VERTS(self.injVerts).V) / self.nbVerts
        expectedPot = ENDTIME * (self.Iinj * self.nbVerts) / (self.capacitance * self.area)

        if MPI.rank == 0:
            print(f'AvgPot = {avgPot*1e3}mV, expectedPot={expectedPot*1e3}mV, difference: {abs(avgPot- expectedPot)*1e3}mV')
            self.assertAlmostEqual(avgPot, expectedPot)

    def test_disttetopsplit_setTriIClamp_n2(self):
        """ Test after ENDTIME that the potential corresponds to what is expected given the current that should be injected"""
        ENDTIME = 0.1

        self.sim.newRun()

        self.sim.membrane.Potential = 0
        self.sim.TRIS(self.mesh.surface).IClamp = self.Iinj

        # Check that we cannot set a current clamp on an internal triangle
        with self.mesh.asLocal():
            innerTri = (self.mesh.tris - self.mesh.surface)[0]
            with self.assertRaises(Exception):
                self.sim.TRI(innerTri).IClamp = self.Iinj

        self.sim.run(ENDTIME)

        avgPot = sum(self.sim.TRIS(self.injTris).V) / self.nbSurfTris
        expectedPot = ENDTIME * (self.Iinj * self.nbSurfTris) / (self.capacitance * self.area)

        if MPI.rank == 0:
            print(f'AvgPot = {avgPot*1e3}mV, expectedPot={expectedPot*1e3}mV, difference: {abs(avgPot- expectedPot)*1e3}mV')
            self.assertAlmostEqual(avgPot, expectedPot)

    def test_disttetopsplit_setVertIClamp_setMembCapac_n2(self):
        """ Test after ENDTIME that the potential corresponds to what is expected given the current that should be injected"""
        ENDTIME = 0.1

        self.sim.newRun()

        self.sim.membrane.Potential = 0
        self.sim.VERTS(self.injVerts).IClamp = self.Iinj
        self.sim.membrane.Capac = 2 * self.capacitance

        self.sim.run(ENDTIME)

        avgPot = sum(self.sim.VERTS(self.injVerts).V) / self.nbVerts
        expectedPot = ENDTIME * (self.Iinj * self.nbVerts) / (2 * self.capacitance * self.area)

        if MPI.rank == 0:
            print(f'AvgPot = {avgPot*1e3}mV, expectedPot={expectedPot*1e3}mV, difference: {abs(avgPot- expectedPot)*1e3}mV')
            self.assertAlmostEqual(avgPot, expectedPot)

    def test_disttetopsplit_setVertIClamp_n1(self):
        # Check that we have the same behavior for 1 rank
        self.test_disttetopsplit_setVertIClamp_n2()

    def test_disttetopsplit_setTriIClamp_n1(self):
        # Check that we have the same behavior for 1 rank
        self.test_disttetopsplit_setTriIClamp_n2()

    def _test_VClamped(self, clampFunc, potFunc, otherPotFunc):
        ENDTIME = 0.1

        self.sim.newRun()

        self.sim.membrane.Potential = 0
        self.sim.VERTS().IClamp = self.Iinj

        clampFunc()

        self.sim.run(ENDTIME)

        avgPotOthers = np.mean(otherPotFunc())
        potClamped = potFunc()

        if MPI.rank == 0:
            self.assertGreater(avgPotOthers, 0)
            self.assertEqual(potClamped, 0)

    def test_disttetopsplit_setVertVClamped_n2(self):
        """ Test after ENDTIME that voltage is clamped in vertex"""
        vert = self.mesh.surface.verts[0]
        def clampFunc():
            self.sim.VERT(vert).VClamped = True
        self._test_VClamped(
            clampFunc,
            lambda:self.sim.VERT(vert).V,
            lambda:self.sim.VERTS(self.mesh.surface.verts - vert.toList()).V
        )

    def test_disttetopsplit_setTriVClamped_n2(self):
        """ Test after ENDTIME that voltage is clamped in triangle"""
        tri = self.mesh.surface[0]
        def clampFunc():
            self.sim.TRI(tri).VClamped = True
        self._test_VClamped(
            clampFunc,
            lambda:self.sim.TRI(tri).V,
            lambda:self.sim.TRIS(self.mesh.surface - tri.toList()).V
        )

    def test_disttetopsplit_setTetVClamped_n2(self):
        """ Test after ENDTIME that voltage is clamped in tetrahedron"""
        tet = self.mesh.tets[0]
        def clampFunc():
            self.sim.TET(tet).VClamped = True
        self._test_VClamped(
            clampFunc,
            lambda:self.sim.TET(tet).V,
            lambda:[self.sim.TETS(self.mesh.tets - tet.toList()).V]
        )


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(DistTetopsplitEFieldClamps))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
