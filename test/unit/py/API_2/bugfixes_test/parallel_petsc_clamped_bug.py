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

"""Bugfix test for PETSc solvers"""

import os
import unittest

import steps.interface
from steps.saving import *
from steps.sim import *

from steps.geom import *
from steps.model import *
from steps.rng import *

FILEDIR = os.path.dirname(os.path.abspath(__file__))
MESHDIR = os.path.join(FILEDIR, "..", "..", "..", "..", "mesh")


class PetscClampedBugfix(unittest.TestCase):
    """Test that unclamping vertices does not trigger a petsc crash"""

    def setUp(self):
        mdl = Model()
        with mdl:
            S1 = Species.Create()
            vsys = VolumeSystem.Create()
            with vsys:
                Diffusion(S1, 0)

        self.mesh = TetMesh.LoadGmsh(os.path.join(MESHDIR, "box.msh"), 1e-6)
        with self.mesh:
            comp = Compartment.Create(self.mesh.tets, vsys)
            patch = Patch.Create(self.mesh.surface, comp, None)
            memb = Membrane.Create([patch])

        rng = RNG()
        self.sim = Simulation("TetVesicle", mdl, self.mesh, rng, MPI.EF_DV_PETSC)

    def testVClampPetsc(self):
        self.sim.newRun()
        cverts = self.mesh.surface[0].verts
        ncverts = self.mesh.surface.verts - cverts
        self.sim.memb.Potential = 0
        self.sim.TRIS(self.mesh.surface).IClamp = 1e-6
        self.sim.VERTS(cverts).VClamped = True
        self.sim.run(10e-6)
        self.assertEqual(self.sim.VERT(cverts[0]).V, 0)
        self.assertGreater(self.sim.VERT(ncverts[0]).V, 0)
        self.sim.VERTS(cverts).VClamped = False
        self.sim.run(20e-6)
        self.assertGreater(self.sim.VERT(cverts[0]).V, 0)
        self.assertGreater(self.sim.VERT(ncverts[0]).V, 0)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(PetscClampedBugfix))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
