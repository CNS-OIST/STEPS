### #################################################################################
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

""" Unit tests for TetVesicle vesicle mobility behavior."""

import os
import unittest

import numpy as np
import steps.interface
from steps.saving import *
from steps.sim import *
from steps.utils import *
from steps.geom import *
from steps.model import *
from steps.rng import *

TEST_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", ".."
)
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class TetVesicleMobilityTestCase(unittest.TestCase):
    """Test vesicle mobility"""

    ENDT = 1
    ves_diam = 40e-9
    ves_dcst = 0.11e-12
    raft_dcst = 1e-12

    def setUp(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            S1, S2 = Species.Create()
            vsys = VolumeSystem.Create()
            vssys = VesicleSurfaceSystem.Create()
            rssys = RaftSurfaceSystem.Create()
            ves1 = Vesicle.Create(self.ves_diam, self.ves_dcst, vssys)
            raft1 = Raft.Create(self.ves_diam, self.raft_dcst, rssys)

            with vsys:
                Diffusion(S1, 0)

            with vssys:
                S1.v > r[1] > None
                r[1].K = 20
                r[1].Immobilization = IMMOBILIZING

                S2.v > r[2] > None
                r[2].K = 20
                r[2].Immobilization = MOBILIZING

            with rssys:
                S1.r > r[1] > None
                r[1].K = 20
                r[1].Immobilization = IMMOBILIZING

                S2.r > r[2] > None
                r[2].K = 20
                r[2].Immobilization = MOBILIZING

        mesh = TetMesh.LoadGmsh(os.path.join(MESH_DIR, "cube.msh"), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets, vsys)
            patch = Patch.Create(mesh.surface, comp, None)

        rng = RNG("mt19937", 512, 987)
        self.sim = Simulation("TetVesicle", mdl, mesh, rng, MPI.EF_NONE)

    # Vesicle diffusion tests #####################################################################

    def testVesDcst(self):
        mesh = self.sim.geom

        # Default case, normal diffusion
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref.Pos = mesh.bbox.center
        self.assertEqual(vesref.Dcst, self.ves_dcst)
        self.sim.run(self.ENDT)
        self.assertNotEqual(vesref.Pos, list(mesh.bbox.center))

        # Single vesicle dcst set to 0
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1", dcst=0)
        vesref.Pos = mesh.bbox.center
        self.assertEqual(vesref.Dcst, 0)
        self.sim.run(self.ENDT)
        self.assertEqual(vesref.Pos, list(mesh.bbox.center))

        # Setting dcst with custom unit
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1", dcst=Parameter(0.1, "nm^2 s^-1"))
        vesref.Pos = mesh.bbox.center
        self.assertAlmostEqual(vesref.Dcst, 0.1e-18)
        self.sim.run(self.ENDT)
        self.assertLess(np.linalg.norm(vesref.Pos - mesh.bbox.center), 1e-9)

        # Setting dcst after vesicle creation
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref.Pos = mesh.bbox.center
        vesref.Dcst = 0
        self.assertEqual(vesref.Dcst, 0)
        self.sim.run(self.ENDT)
        self.assertEqual(vesref.Pos, list(mesh.bbox.center))

        tet = mesh.tets[mesh.bbox.center]
        # Tet dcst set to 0
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref.Pos = mesh.bbox.center
        self.assertEqual(self.sim.TET(tet).ves1.Dcst, -1)
        self.sim.TET(tet).ves1.Dcst = 0
        self.assertEqual(self.sim.TET(tet).ves1.Dcst, 0)
        self.assertFalse(self.sim.TET(tet).ves1.DcstRel)
        self.sim.run(self.ENDT)
        self.assertEqual(vesref.Pos, list(mesh.bbox.center))

        # Vesicles dcst = 0 and Tet dcst > 0
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1", dcst=0)
        vesref.Pos = mesh.bbox.center
        self.assertEqual(self.sim.TET(tet).ves1.Dcst, -1)
        self.sim.TET(tet).ves1.Dcst = self.ves_dcst
        self.assertEqual(self.sim.TET(tet).ves1.Dcst, self.ves_dcst)
        self.assertFalse(self.sim.TET(tet).ves1.DcstRel)
        self.sim.run(self.ENDT)
        self.assertNotEqual(vesref.Pos, list(mesh.bbox.center))

        # Tet dcst set relative to vesicle dcst
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref.Pos = mesh.bbox.center
        self.assertEqual(self.sim.TET(tet).ves1.Dcst, -1)
        self.sim.TET(tet).ves1.Dcst = Params(1e-50, relative=True)
        self.assertEqual(self.sim.TET(tet).ves1.Dcst, 1e-50)
        self.assertTrue(self.sim.TET(tet).ves1.DcstRel)
        self.sim.run(self.ENDT)
        self.assertAlmostEqual(np.linalg.norm(vesref.Pos - mesh.bbox.center), 0)

    # Vesicle mobility tests ######################################################################

    def testVesGetSetMobility(self):
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        self.assertEqual(vesref.Immobility, 0)

        vesref.Immobility = 1
        self.assertEqual(vesref.Immobility, 1)

        vesref.Immobility = 10
        self.assertEqual(vesref.Immobility, 10)

        vesref.Immobility = 0
        self.assertEqual(vesref.Immobility, 0)

        with self.assertRaises(Exception):
            vesref.Immobility = -1

        with self.assertRaises(Exception):
            vesref.Immobility = "1"

    def testVesMobility(self):
        # Unrestrained
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        initPos = vesref.Pos
        self.sim.run(self.ENDT)

        self.assertNotEqual(vesref.Pos, initPos)

        # Restrained
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        initPos = vesref.Pos
        vesref.Immobility = 1
        self.sim.run(self.ENDT)

        self.assertEqual(vesref.Pos, initPos)

    def testVesMobilization(self):
        # Immobilizing reaction
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref("surf").S1.Count = 5
        initPos = vesref.Pos
        self.sim.run(self.ENDT)
        self.assertGreater(vesref.Immobility, 0)
        self.assertNotEqual(vesref.Pos, initPos)

        # Mobilizing reaction starting unrestrained
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref("surf").S2.Count = 5
        initPos = vesref.Pos
        self.sim.run(self.ENDT)
        self.assertEqual(vesref.Immobility, 0)
        self.assertNotEqual(vesref.Pos, initPos)

        # Mobilizing reaction starting restrained
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref("surf").S2.Count = 5
        vesref.Immobility = 5
        initPos = vesref.Pos
        self.sim.run(self.ENDT)
        self.assertEqual(vesref.Immobility, 0)
        self.assertNotEqual(vesref.Pos, initPos)

        # Mobilizing reaction starting restrained 2
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref("surf").S2.Count = 5
        vesref.Immobility = 6
        initPos = vesref.Pos
        self.sim.run(self.ENDT)
        self.assertEqual(vesref.Immobility, 1)
        self.assertEqual(vesref.Pos, initPos)

        # Both mobilizing and immobilizing reactions starting unrestrained
        self.sim.newRun()
        vesref = self.sim.comp.addVesicle("ves1")
        vesref("surf").S1.Count = 5
        vesref("surf").S2.Count = 5
        initPos = vesref.Pos
        self.sim.run(self.ENDT)
        self.assertGreaterEqual(vesref.Immobility, 0)
        self.assertNotEqual(vesref.Pos, initPos)

    # Raft mobility tests #########################################################################

    def testRaftGetSetMobility(self):
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        self.assertEqual(raftref.Immobility, 0)

        raftref.Immobility = 1
        self.assertEqual(raftref.Immobility, 1)

        raftref.Immobility = 10
        self.assertEqual(raftref.Immobility, 10)

        raftref.Immobility = 0
        self.assertEqual(raftref.Immobility, 0)

        with self.assertRaises(Exception):
            raftref.Immobility = -1

        with self.assertRaises(Exception):
            raftref.Immobility = "1"

    def testRaftMobility(self):
        # Unrestrained
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        initPos = raftref.Pos
        self.sim.run(self.ENDT)

        self.assertNotEqual(raftref.Pos, initPos)

        # Restrained
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        initPos = raftref.Pos
        raftref.Immobility = 1
        self.sim.run(self.ENDT)

        self.assertEqual(raftref.Pos, initPos)

    def testRaftMobilization(self):
        # Immobilizing reaction
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        raftref.S1.Count = 5
        initPos = raftref.Pos
        self.sim.run(self.ENDT)
        self.assertGreater(raftref.Immobility, 0)
        self.assertNotEqual(raftref.Pos, initPos)

        # Mobilizing reaction starting unrestrained
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        raftref.S2.Count = 5
        initPos = raftref.Pos
        self.sim.run(self.ENDT)
        self.assertEqual(raftref.Immobility, 0)
        self.assertNotEqual(raftref.Pos, initPos)

        # Mobilizing reaction starting restrained
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        raftref.S2.Count = 5
        raftref.Immobility = 5
        initPos = raftref.Pos
        self.sim.run(self.ENDT)
        self.assertEqual(raftref.Immobility, 0)
        self.assertNotEqual(raftref.Pos, initPos)

        # Mobilizing reaction starting restrained 2
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        raftref.S2.Count = 5
        raftref.Immobility = 6
        initPos = raftref.Pos
        self.sim.run(self.ENDT)
        self.assertEqual(raftref.Immobility, 1)
        self.assertEqual(raftref.Pos, initPos)

        # Both mobilizing and immobilizing reactions starting unrestrained
        self.sim.newRun()
        raftref = self.sim.TRI(self.sim.patch.tris[0]).addRaft("raft1")
        raftref.S1.Count = 5
        raftref.S2.Count = 5
        initPos = raftref.Pos
        self.sim.run(self.ENDT)
        self.assertGreaterEqual(raftref.Immobility, 0)
        self.assertNotEqual(raftref.Pos, initPos)


def suite():
    all_tests = []
    all_tests.append(
        unittest.TestLoader().loadTestsFromTestCase(TetVesicleMobilityTestCase)
    )
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
