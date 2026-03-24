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

""" Unit tests for TetVesicle LinkSpecies constraints."""

import os
import unittest
import numpy as np

from steps import interface
from steps.saving import *
from steps.sim import *
from steps.geom import *
from steps.model import *
from steps.rng import *

TEST_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "..", ".."
)
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class TetVesicleVesBindingTestCase(unittest.TestCase):
    """Test vesicle paths"""

    ENDT = 1

    nb_vesicles = 200
    ves_diam = 40e-9
    ves_dcst = 0.11e-12

    ls_dcst = 0.1e-12

    vesbind_kcst = 1e9
    min_length = 10e-9
    max_length = 30e-9
    max_angle_L1 = np.pi
    max_angle_L2 = np.pi / 8

    min_link_number = 20

    def setUp(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            S1, S2 = Species.Create()
            L1 = LinkSpecies.Create(dcst=self.ls_dcst)
            L2 = LinkSpecies.Create(dcst=self.ls_dcst, max_angle=self.max_angle_L2)

            self.assertEqual(L1.MaxAngle, self.max_angle_L1)
            self.assertEqual(L2.MaxAngle, self.max_angle_L2)

            vsys = VolumeSystem.Create()
            ves1 = Vesicle.Create(self.ves_diam, self.ves_dcst)

            with vsys:
                ves1bind1 = VesicleBind.Create(
                    (ves1, ves1),
                    (S1, S1),
                    (L1, L1),
                    self.min_length,
                    self.max_length,
                    kcst=self.vesbind_kcst,
                )
                ves1bind2 = VesicleBind.Create(
                    (ves1, ves1),
                    (S2, S2),
                    (L2, L2),
                    self.min_length,
                    self.max_length,
                    kcst=self.vesbind_kcst,
                )

                ves1bind3 = VesicleBind.Create(
                    (ves1, ves1),
                    (S1, S2),
                    (L1, L2),
                    self.min_length,
                    self.max_length,
                    kcst=self.vesbind_kcst,
                )

        mesh = TetMesh.LoadGmsh(os.path.join(MESH_DIR, "cube.msh"), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets, vsys)

        rng = RNG("mt19937", 512, 987)
        self.sim = Simulation("TetVesicle", mdl, mesh, rng, MPI.EF_NONE)

    def testVesicleBindingConstraints(self):
        sim = self.sim
        mesh = sim.geom
        mdl = sim.model

        sim.newRun()

        sim.comp.ves1.Count = self.nb_vesicles
        sim.comp.VESICLES(mdl.ves1)("surf").ALL(Species).Count = 10

        sim.run(self.ENDT)

        ves_pos = sim.comp.VESICLES(mdl.ves1).Pos

        L1_pos = sim.comp.VESICLES(mdl.ves1)("surf").LINKSPECS(mdl.L1).Pos
        L1_links = sim.comp.VESICLES(mdl.ves1)("surf").LINKSPECS(mdl.L1).LinkedTo

        L2_pos = sim.comp.VESICLES(mdl.ves1)("surf").LINKSPECS(mdl.L2).Pos
        L2_links = sim.comp.VESICLES(mdl.ves1)("surf").LINKSPECS(mdl.L2).LinkedTo

        # Build link spec dict
        L1_ves = {
            ls_ind: ves_ind for ves_ind, ls_dct in L1_pos.items() for ls_ind in ls_dct
        }
        L2_ves = {
            ls_ind: ves_ind for ves_ind, ls_dct in L2_pos.items() for ls_ind in ls_dct
        }

        # First test L1 links, should default to max angle being pi
        angles = []
        for ves_ind, ls_dct in L1_pos.items():
            vpos = np.array(ves_pos[ves_ind])
            for l1, p1 in ls_dct.items():
                l2 = L1_links[ves_ind][l1]
                p1 = np.array(p1)
                if l2 in L1_ves:
                    p2 = np.array(L1_pos[L1_ves[l2]][l2])
                else:
                    p2 = np.array(L2_pos[L2_ves[l2]][l2])
                A = p1 - vpos
                B = p2 - p1
                dist = np.linalg.norm(B)
                angle = np.arccos(A.dot(B) / (np.linalg.norm(A) * dist))
                angles.append(angle)

                self.assertGreaterEqual(dist, self.min_length)
                self.assertLessEqual(dist, self.max_length)
                self.assertLessEqual(angle, self.max_angle_L1)
        # Check that some links were formed
        self.assertGreater(len(angles), self.min_link_number)

        # Then test L2 links
        angles = []
        for ves_ind, ls_dct in L2_pos.items():
            vpos = np.array(ves_pos[ves_ind])
            for l1, p1 in ls_dct.items():
                l2 = L2_links[ves_ind][l1]
                p1 = np.array(p1)
                if l2 in L1_ves:
                    p2 = np.array(L1_pos[L1_ves[l2]][l2])
                else:
                    p2 = np.array(L2_pos[L2_ves[l2]][l2])
                A = p1 - vpos
                B = p2 - p1
                dist = np.linalg.norm(B)
                angle = np.arccos(A.dot(B) / (np.linalg.norm(A) * dist))
                angles.append(angle)

                self.assertGreaterEqual(dist, self.min_length)
                self.assertLessEqual(dist, self.max_length)
                self.assertLessEqual(angle, self.max_angle_L2)
        # Check that some links were formed
        self.assertGreater(len(angles), self.min_link_number)


def suite():
    all_tests = []
    all_tests.append(
        unittest.TestLoader().loadTestsFromTestCase(TetVesicleVesBindingTestCase)
    )
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
