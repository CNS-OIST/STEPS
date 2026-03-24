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

""" Unit tests for DistTetOpSplit reaction extent querrying"""

import importlib
import os
import unittest

from steps import interface

from steps.model import *
from steps.sim import *
from steps.geom import *
from steps.rng import *

TEST_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", ".."
)
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class DistTetOpSplitReactionExtent(unittest.TestCase):
    """Test reaction extent querrying of DistTetOpSplit"""

    def setUp(self):
        self.initCount = 200
        self.K = 50
        self.vdK = VDepRate(lambda V: self.K)
        self.ENDT = 1e-3

        mdl = Model()
        r = ReactionManager()

        with mdl:
            vsys = VolumeSystem()
            ssys = SurfaceSystem()
            SA, SB, SC, SD = Species.Create()
            susA, susB, susC, susD = SubUnitState.Create()
            cmplx = Complex.Create([susA, susB, susC, susD], statesAsSpecies=False)

            with vsys:
                # Volume reaction
                SA > r["reac"] > SB
                r["reac"].K = self.K

                # Volume complex reaction
                with cmplx[...]:
                    susA > r["creac"] > susB
                    r["creac"].K = self.K

            with ssys:
                # Surface reaction
                SA.s > r["sreac"] > SB.s
                r["sreac"].K = self.K

                # Complex surface reaction
                with cmplx[...]:
                    susA.s > r["csreac"] > susB.s
                    r["csreac"].K = self.K

                # Voltage-dependent Surface reaction
                SC.s > r["vsreac"] > SD.s
                r["vsreac"].K = self.vdK

                # Voltage-dependent Complex surface reaction
                with cmplx[...]:
                    susC.s > r["vcsreac"] > susD.s
                    r["vcsreac"].K = self.vdK

        mesh = DistMesh(os.path.join(MESH_DIR, "box.msh"), 1e-6)
        with mesh.asLocal():
            comp = Compartment.Create(mesh.tets, vsys, conductivity=1)
            patch = Patch.Create(mesh.surface, comp, None, ssys)
            memb = Membrane.Create([patch], capacitance=1)

        rng = RNG()
        self.sim = Simulation("DistTetOpSplit", mdl, mesh, rng)

    def testReactionExtent_n1(self):
        susA = self.sim.model.susA
        susC = self.sim.model.susC
        paths = [
            (self.sim.comp.reac, self.sim.comp.SB),
            (self.sim.comp.creac, self.sim.comp.cmplx.susB),
            (self.sim.patch.sreac, self.sim.patch.SB),
            (self.sim.patch.csreac, self.sim.patch.cmplx.susB),
            (self.sim.patch.vsreac, self.sim.patch.SD),
            (self.sim.patch.vcsreac, self.sim.patch.cmplx.susD),
        ]

        for ext, cnt in paths:
            with self.subTest(ext=ext, cnt=cnt):
                self.assertEqual(ext.Extent, 0)

        self.sim.newRun()

        self.sim.comp.SA.Count = self.initCount
        self.sim.comp.cmplx[susA].Count = self.initCount
        self.sim.patch.SA.Count = self.initCount
        self.sim.patch.SC.Count = self.initCount
        self.sim.patch.cmplx[susA].Count = self.initCount
        self.sim.patch.cmplx[susC].Count = self.initCount

        self.sim.run(self.ENDT)

        for ext, cnt in paths:
            with self.subTest(ext=ext, cnt=cnt):
                self.assertEqual(ext.Extent, cnt.Count)

        self.sim.newRun()

        for ext, cnt in paths:
            with self.subTest(ext=ext, cnt=cnt):
                self.assertEqual(ext.Extent, 0)


def suite():
    all_tests = []
    all_tests.append(
        unittest.TestLoader().loadTestsFromTestCase(DistTetOpSplitReactionExtent)
    )
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
