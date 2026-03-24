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

"""Bugfixes test for new compartment clamping behaviour"""

import os
import unittest

import numpy as np
import steps.interface
from steps.saving import ResultSelector
from steps.sim import Simulation, SSAMethod, NextEventSearchMethod, DiffusionMethod

from steps.geom import Compartment, DiffBoundary, DistMesh, TetList, TetMesh
from steps.model import Diffusion, Model, ReactionManager, Species, VolumeSystem
from steps.rng import RNG

FILEDIR = os.path.dirname(os.path.abspath(__file__))
MESHDIR = os.path.join(FILEDIR, "..", "..", "..", "..", "mesh")

DT = 0.01
ENDT = 0.1
DCST = 1e-12
INJ = 1000


class CompartmentClamping(unittest.TestCase):
    def getSim(self, solver):
        r = ReactionManager()
        mdl = Model()
        with mdl:
            SA, SB = Species.Create()
            vsys = VolumeSystem.Create()
            with vsys:
                SB > r["r1"] > SA
                r["r1"].K = 0

                Diffusion(SA, DCST)

        if solver.startswith("DistTetOpSplit"):
            mesh = DistMesh(os.path.join(MESHDIR, "cube.msh"), 1e-6)
        else:
            mesh = TetMesh.LoadGmsh(os.path.join(MESHDIR, "cube.msh"), 1e-6)

        with mesh:
            lefttets = TetList(
                [tet for tet in mesh.tets if tet.center.x < mesh.bbox.center.x]
            )
            comp1 = Compartment.Create(lefttets, vsys)
            comp2 = Compartment.Create(mesh.tets - lefttets, vsys)

            diffb = DiffBoundary.Create(comp1.surface & comp2.surface)

        rng = RNG(seed=1234)
        solvArgs = {}
        if solver == "DistTetOpSplitDiffLeap":
            solvArgs["SSAMethod"] = SSAMethod.RLEAPING
            solvArgs["searchMethod"] = NextEventSearchMethod.RLEAPING
            solvArgs["diffMethod"] = DiffusionMethod.TAU_LEAPING_DT
            solver = "DistTetOpSplit"
        sim = Simulation(solver, mdl, mesh, rng, **solvArgs)
        return sim

    def testCompartmentClamping_n4(self):
        solvers = ["Tetexact", "TetOpSplit", "TetVesicle", "DistTetOpSplit", "DistTetOpSplitDiffLeap"]
        for solver in solvers:
            with self.subTest(solver=solver):
                sim = self.getSim(solver)

                rs = ResultSelector(sim)
                counts1 = rs.TETS(sim.geom.comp1.tets).SA.Count
                count2 = rs.comp2.SA.Count
                count1SB = rs.comp1.SB.Count

                sim.toSave(counts1, count2, count1SB, dt=DT)

                # Unclamped compartment, inactive diffusion boundary
                sim.newRun()
                sim.comp1.SA.Count = INJ
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] != counts1.data[-1, -1, :]).any()
                )
                self.assertEqual(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertEqual(count2.data[-1, -1, 0], 0)

                # Unclamped compartment, active diffusion boundary
                sim.newRun()
                sim.diffb.SA.DiffusionActive = True
                sim.comp1.SA.Count = INJ
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] != counts1.data[-1, -1, :]).any()
                )
                self.assertLess(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertGreater(count2.data[-1, -1, 0], 0)
                self.assertEqual(
                    np.sum(counts1.data[-1, -1, :]) + count2.data[-1, -1, 0], INJ
                )

                # Clamped compartment, inactive diffusion boundary
                sim.newRun()
                sim.comp1.SA.Clamped = True
                sim.comp1.SA.Count = INJ
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] != counts1.data[-1, -1, :]).any()
                )
                self.assertEqual(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertEqual(count2.data[-1, -1, 0], 0)

                # Clamped compartment, active diffusion boundary
                sim.newRun()
                sim.comp1.SA.Clamped = True
                sim.diffb.SA.DiffusionActive = True
                sim.comp1.SA.Count = INJ
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] != counts1.data[-1, -1, :]).any()
                )
                self.assertEqual(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertGreater(count2.data[-1, -1, 0], 0)
                self.assertGreater(
                    np.sum(counts1.data[-1, -1, :]) + count2.data[-1, -1, 0], INJ
                )

                # Clamped tets, inactive diffusion boundary
                sim.newRun()
                sim.TETS(sim.geom.comp1.tets).SA.Clamped = True
                sim.comp1.SA.Count = INJ
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] == counts1.data[-1, -1, :]).all()
                )
                self.assertEqual(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertEqual(count2.data[-1, -1, 0], 0)

                # Clamped tets, active diffusion boundary
                sim.newRun()
                sim.TETS(sim.geom.comp1.tets).SA.Clamped = True
                sim.diffb.SA.DiffusionActive = True
                sim.comp1.SA.Count = INJ
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] == counts1.data[-1, -1, :]).all()
                )
                self.assertEqual(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertGreater(count2.data[-1, -1, 0], 0)
                self.assertGreater(
                    np.sum(counts1.data[-1, -1, :]) + count2.data[-1, -1, 0], INJ
                )

                # Clamped compartment, activated reaction
                sim.newRun()
                sim.comp1.SA.Clamped = True
                sim.comp1.SA.Count = INJ
                sim.comp1.SB.Count = INJ
                sim.comp1.r1.K = 100
                sim.run(ENDT)
                self.assertTrue(
                    (counts1.data[-1, 0, :] != counts1.data[-1, -1, :]).any()
                )
                self.assertEqual(np.sum(counts1.data[-1, -1, :]), INJ)
                self.assertEqual(count2.data[-1, -1, 0], 0)
                self.assertLess(count1SB.data[-1, -1, 0], INJ)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(CompartmentClamping))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
