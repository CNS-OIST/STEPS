### #################################################################################
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

""" Unit tests for TetVesicle vesicle paths behavior."""

import numpy as np
import os
import unittest

from steps import interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.sim import *
from steps.saving import *


TEST_DIR = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), "..", "..", "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")


class TetVesiclePathsTestCase(unittest.TestCase):
    """Test vesicle paths"""

    ENDT = 1
    ves_diam = 40e-9
    ves_dcst = 0.11e-12

    def setUp(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            S1 = Species.Create()
            vsys = VolumeSystem.Create()
            ves1 = Vesicle.Create(self.ves_diam, self.ves_dcst)

            with vsys:
                Diffusion(S1, 0)

        mesh = TetMesh.LoadGmsh(os.path.join(MESH_DIR, 'cube.msh'), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets, vsys)

        rng = RNG('mt19937', 512, 987)
        self.sim = Simulation('TetVesicle', mdl, mesh, rng, MPI.EF_NONE)

    def addPath(self, speed, step_size):
        path = self.sim.addVesiclePath(f'path')
        bbox = self.sim.geom.bbox
        p1 = path.addPoint(bbox.center)
        p2 = path.addPoint([bbox.max.x, bbox.center.y, bbox.center.z])
        path.addBranch(p1, {p2: 1})
        path.addVesicle(self.sim.model.ves1, speed=speed,
                        stoch_stepsize=step_size)

    def testPathMeshBoundaryCentered(self):
        step_size = 1e-9
        self.addPath(3e-6, step_size)
        sim = self.sim
        mesh = sim.geom

        sim.newRun()
        sim.setVesicleDT(1e-2)
        vesref = sim.comp.addVesicle('ves1')
        vesref.Pos = mesh.bbox.center

        sim.run(self.ENDT)
        dist_to_boundary = (
            mesh.bbox.max.x - self.ves_diam / 2) - vesref.Pos[0]
        self.assertLess(dist_to_boundary, step_size)

    def testPathMeshBoundaryOffCentered(self):
        step_size = 1e-9
        self.addPath(3e-6, step_size)
        sim = self.sim
        mesh = sim.geom

        sim.newRun()
        sim.setVesicleDT(1e-2)
        vesref = sim.comp.addVesicle('ves1')
        vesref.Pos = mesh.bbox.center + np.array([5.12e-9, 3.45e-9, 2.78e-9])

        sim.run(self.ENDT)
        dist_to_boundary = (
            mesh.bbox.max.x - self.ves_diam / 2) - vesref.Pos[0]
        self.assertLess(dist_to_boundary, step_size)

    def testVesicleSpeed(self):
        speed = 1e-7
        self.addPath(speed, 1e-9)
        speed_tolerance = 5e-2

        sim = self.sim
        mesh = sim.geom

        sim.newRun()
        vesref = sim.comp.addVesicle('ves1')
        vesref.Pos = mesh.bbox.center

        end_t = (mesh.bbox.max.x - mesh.bbox.center.x -
                 self.ves_diam / 2) / speed * 0.1
        sim.run(end_t)
        exp_speed = (vesref.Pos[0] - mesh.bbox.center.x) / end_t
        self.assertLess(
            abs(exp_speed - speed) / speed, speed_tolerance)

    def testPathEnd(self):
        path = self.sim.addVesiclePath(f'path')
        bbox = self.sim.geom.bbox
        p1 = path.addPoint(bbox.center)
        shift = (bbox.max.x - bbox.min.x) / 100
        nbPoints = 10
        for i in range(nbPoints):
            p2 = path.addPoint([bbox.center.x + (i + 1) / nbPoints * shift, bbox.center.y, bbox.center.z])
            path.addBranch(p1, {p2: 1})
            p1 = p2
        path.addVesicle(self.sim.model.ves1, speed=1e-7, stoch_stepsize=1e-10)

        sim = self.sim
        mesh = sim.geom

        rs = ResultSelector(sim)
        pos = rs.VESICLES(sim.comp.ves1).Pos
        sim.toSave(pos, dt=0.01)

        sim.newRun()
        vesref = sim.comp.addVesicle('ves1')
        vesref.Pos = mesh.bbox.center

        sim.run(self.ENDT)

        self.assertGreater(max([np.linalg.norm(p - mesh.bbox.center) for dct in pos.data[0, :, 0] for idx, p in dct.items()]), 2 * shift)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetVesiclePathsTestCase))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
