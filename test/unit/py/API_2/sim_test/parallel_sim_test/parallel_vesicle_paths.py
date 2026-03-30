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
            ves2 = Vesicle.Create(self.ves_diam, 0)

            with vsys:
                Diffusion(S1, 0)

        mesh = TetMesh.LoadGmsh(os.path.join(MESH_DIR, 'cube.msh'), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets, vsys)

        rng = RNG('mt19937', 512, 987)
        self.sim = Simulation('TetVesicle', mdl, mesh, rng, MPI.EF_NONE)

    def addPath(self, speed, step_size, ves_type='ves1', bind_to_start=True, allow_edge_binding=False, **kwargs):
        path = self.sim.addVesiclePath(f'path', bind_to_start=bind_to_start)
        bbox = self.sim.geom.bbox
        p1 = path.addPoint(bbox.center)
        p2 = path.addPoint([bbox.max.x, bbox.center.y, bbox.center.z])
        path.addEdge(p1, p2, allow_binding=allow_edge_binding)
        path.addVesicle(ves_type, speed=speed, stoch_stepsize=step_size, **kwargs)

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
        self.addPath(speed, 1e-10)
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
            path.addEdge(p1, p2)
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

    def _testPathEdgeBinding(self, rate=1.0, min_br=0, max_br=None, ves_type='ves2', ves_pos=None, expected_rate=0):
        """Test that the observed rate of binding to a path edges matches what is expected"""
        sim = self.sim
        mesh = sim.geom

        # Create a path across the mesh along the x axis
        path = sim.addVesiclePath('path', bind_to_start=False)
        c = mesh.bbox.center
        p1 = path.addPoint([mesh.bbox.min.x, c.y, c.z])
        p2 = path.addPoint([mesh.bbox.max.x, c.y, c.z])
        path.addEdge(p1, p2, allow_binding=True)
        path.addVesicle(
            ves_type, speed=1e-7,
            binding_rate=rate, min_binding_radius=min_br, max_binding_radius=max_br
        )

        # Run the simulation and regularly check whether the vesicle is bound to the path
        # If so, remove it and add another vesicle, keeping track of the time between bind events
        sim.newRun()
        vesref = sim.comp.addVesicle(ves_type)
        if ves_pos is not None:
            vesref.Pos = ves_pos
        last_binding = 0
        binding_times = []
        if hasattr(expected_rate, '__call__'):
            expected_rate = expected_rate(mesh)
        if expected_rate > 0:
            ENDT = 1000 / expected_rate
            DT = 0.1 / expected_rate
        else:
            ENDT = 1
            DT = 1e-3
        for t in np.arange(DT, ENDT, DT):
            sim.run(t)
            if vesref.OnPath is not None:
                binding_times.append(t - last_binding)
                last_binding = t

                self.assertEqual(vesref.OnPath[0], 'path')

                vesref.delete()
                vesref = sim.comp.addVesicle(ves_type)
                if ves_pos is not None:
                    vesref.Pos = ves_pos

        if len(binding_times) > 0:
            measured_rate = 1 / np.mean(binding_times)
        else:
            measured_rate = 0

        if MPI.rank == 0:
            print(f'{measured_rate=}, {expected_rate=}')
        self.assertTrue(np.isclose(measured_rate, expected_rate, rtol=0.1, atol=0))

    def testPathEdgeBinding(self):
        mesh = self.sim.geom
        params = {
            'Default values, centered':
                dict(
                    rate=1e3, min_br=0, max_br=-1,
                    ves_pos=mesh.bbox.center,
                    expected_rate=1e3,
                ),
            'Default values, off-centered, further than binding radius':
                dict(
                    rate=1e3, min_br=0, max_br=-1,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam / 2 + 1e-10, 0]),
                    expected_rate=0,
                ),
            'Default values, off-centered, within binding radius':
                dict(
                    rate=1e3, min_br=0, max_br=-1,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam / 4, 0]),
                    expected_rate=1e3,
                ),

            'Higher maximum binding radius, centered':
                dict(
                    rate=1e3, min_br=0, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center,
                    expected_rate=1e3,
                ),
            'Higher maximum binding radius, off-centered, further than binding radius':
                dict(
                    rate=1e3, min_br=0, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam + 1e-10, 0]),
                    expected_rate=0,
                ),
            'Higher maximum binding radius, off-centered, within binding radius':
                dict(
                    rate=1e3, min_br=0, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam / 2, 0]),
                    expected_rate=1e3,
                ),

            'Higher binding radius, no binding inside vesicle, centered':
                dict(
                    rate=1e3, min_br=self.ves_diam / 2, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center,
                    expected_rate=0,
                ),
            'Higher binding radius, no binding inside vesicle, off-centered, further than binding radius':
                dict(
                    rate=1e3, min_br=self.ves_diam / 2, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam + 1e-10, 0]),
                    expected_rate=0,
                ),
            'Higher binding radius, no binding inside vesicle, off-centered, within forbidden binding radius':
                dict(
                    rate=1e3, min_br=self.ves_diam / 2, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam / 3, 0]),
                    expected_rate=0,
                ),
            'Higher binding radius, no binding inside vesicle, off-centered, within binding radius':
                dict(
                    rate=1e3, min_br=self.ves_diam / 2, max_br=self.ves_diam,
                    ves_pos=mesh.bbox.center + np.array([0, self.ves_diam * 3 / 4, 0]),
                    expected_rate=1e3,
                ),
            'Low binding rate, centered':
                dict(
                    rate=1e1, min_br=0, max_br=-1,
                    ves_pos=mesh.bbox.center,
                    expected_rate=1e1,
                ),
        }
        for name, kwargs in params.items():
            # Try with absolute values
            self.setUp()
            with self.subTest(name, **kwargs):
                self._testPathEdgeBinding(**kwargs)
            self.tearDown()

            # Try with relative values
            self.setUp()
            if kwargs['min_br'] > 0:
                kwargs['min_br'] = - kwargs['min_br'] / (self.ves_diam / 2)
            if kwargs['max_br'] > 0:
                kwargs['max_br'] = - kwargs['max_br'] / (self.ves_diam / 2)
            with self.subTest(name, **kwargs):
                self._testPathEdgeBinding(**kwargs)
            self.tearDown()

    def testPathEdgeUnbinding(self, expected_rate=50):
        self.addPath(
            1e-10, 1e-9, ves_type='ves2', bind_to_start=False, allow_edge_binding=True,
            binding_rate=1e6, unbinding_rate=expected_rate
        )

        sim = self.sim
        mesh = sim.geom

        sim.newRun()
        vesref = sim.comp.addVesicle('ves2')
        vesref.Pos = mesh.bbox.center

        DT = 1e-3
        ENDT = 10
        bound = False
        last_binding = None
        unbinding_times = []
        for t in np.arange(DT, ENDT, DT):
            sim.run(t)
            onPath = vesref.OnPath
            if not bound and onPath is not None:
                bound = True
                last_binding = t
                self.assertEqual(onPath[0], 'path')
            elif bound and onPath is None:
                bound = False
                unbinding_times.append(t - last_binding)

        if len(unbinding_times) > 0:
            measured_rate = 1 / np.mean(unbinding_times)
        else:
            measured_rate = 0

        if MPI.rank == 0:
            print(f'{measured_rate=}, {expected_rate=}')
        self.assertTrue(np.isclose(measured_rate, expected_rate, rtol=0.1, atol=0))

    def _testPathIntersection(self, allow_path_intersection, expected_pos, abs_tol):
        self.addPath(
            1e-7, 1e-10, ves_type='ves2', bind_to_start=True,
            min_binding_radius=self.ves_diam/2, max_binding_radius=self.ves_diam,
            allow_path_intersection=allow_path_intersection
        )

        sim = self.sim
        mesh = sim.geom

        sim.newRun()
        vesref = sim.comp.addVesicle('ves2')
        vesref.Pos = mesh.bbox.center - np.array([self.ves_diam / 2 + 1e-9, 0, 0])

        sim.run(self.ENDT)
        self.assertEqual(vesref.OnPath[0], 'path')
        pos = vesref.Pos
        if MPI.rank == 0:
            print(pos, f'{expected_pos=}')
        self.assertLessEqual(np.linalg.norm(pos - expected_pos), abs_tol)

    def testPathIntersection(self):
        with self.subTest('Allow intersection'):
            mesh = self.sim.geom
            self._testPathIntersection(True, mesh.bbox.center + np.array([self.ENDT * 1e-7 - self.ves_diam / 2 - 1e-9, 0, 0]), 2e-9)
        self.tearDown()
        self.setUp()
        with self.subTest('Prevent intersection'):
            mesh = self.sim.geom
            self._testPathIntersection(False, mesh.bbox.center - np.array([self.ves_diam / 2, 0, 0]), 1e-10)



def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetVesiclePathsTestCase))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
