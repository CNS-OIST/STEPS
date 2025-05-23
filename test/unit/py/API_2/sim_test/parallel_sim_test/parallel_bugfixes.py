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

"""Small tests for various bugfixes"""

import os
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

MESHDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../../mesh')

class TestTetOpSplitOwnedClamped(unittest.TestCase):

    def testSetClamped(self):
        mdl = Model()
        with mdl:
            S1 = Species.Create()
            vsys = VolumeSystem()
            ssys = SurfaceSystem()
            with vsys:
                Diffusion(S1, 1e-12)
            with ssys:
                Diffusion(S1, 1e-12)

        mesh = TetMesh.LoadGmsh(os.path.join(MESHDIR, 'cube.msh'), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets, vsys)
            patch = Patch.Create(mesh.surface, comp, None, ssys)

        rng = RNG('mt19937', 1024, 1234)
        part = MeshPartition([0, 1])
        sim = Simulation('TetOpSplit', mdl, mesh, rng, False)

        sim.newRun()

        sim.LIST(comp, patch).S1.Count = 10
        sim.LIST(comp, patch).S1.Clamped = True

        sim.run(1)

        self.assertEqual(sim.comp.S1.Count, 10)
        self.assertEqual(sim.patch.S1.Count, 10)


class VesicleDiffusionCompartments(unittest.TestCase):
    def testVesicleDiffForbiddenComp(self):
        DCST = 1e-13
        VES_DIAM = 40e-9
        VES_NB = 50
        SAVE_DT = 0.001

        mdl = Model()
        with mdl:
            S1 = Species.Create()
            vsys = VolumeSystem()

            ves1 = Vesicle.Create(VES_DIAM, DCST)
            with vsys:
                Diffusion(S1, 0)

        mesh = TetMesh.LoadGmsh(os.path.join(MESHDIR, 'box.msh'), 1e-6)
        with mesh:
            comp1 = Compartment.Create(
                [tet for tet in mesh.tets if tet.center.x > mesh.bbox.center.x], vsys)
            comp2 = Compartment.Create(mesh.tets - comp1.tets, vsys)

        rng = RNG()
        sim = Simulation('TetVesicle', mdl, mesh, rng, MPI.EF_NONE)

        sim.newRun()
        sim.comp1.ves1.Count = VES_NB
        sim.run(1)
        self.assertEqual(sim.comp1.ves1.Count, VES_NB)
        self.assertEqual(sim.comp2.ves1.Count, 0)

