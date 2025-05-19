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

""" Unit tests for vesicle geometry and related methods."""

import unittest
import tempfile
import os

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.utils import *

FILEDIR = os.path.dirname(os.path.abspath(__file__))

class vesicleGeom(unittest.TestCase):
    """Test tetmesh loading and geometrical elements handling."""
    def setUp(self):
        self.mesh = TetMesh.LoadAbaqus(os.path.join(FILEDIR, 'meshes', 'brick_40_4_4_1400tets.inp'), 1e-6)

    def testEndocyticZoneCreation(self):
        with self.mesh:
            comp1 = Compartment.Create(self.mesh.tets)

            surf = self.mesh.surface

            patch1 = Patch.Create(surf[:len(surf)//2], comp1, None)

            with patch1:
                triLst = patch1.tris[:len(patch1.tris)//2]
                endoz1 = EndocyticZone.Create(triLst)

            # Without using a patch
            with self.assertRaises(Exception):
                endoz2 = EndocyticZone.Create(triLst)

            with self.assertRaises(TypeError):
                with patch1:
                    endosz3 = EndocyticZone.Create(42)

            with self.assertRaises(Exception):
                with patch1:
                    endosz4 = EndocyticZone.Create(surf[len(surf)//2:])

            # From steps object
            so = endoz1._getStepsObjects()[0]
            so1 = so.__class__('endoz1b', patch1._getStepsObjects()[0], triLst.indices)
            with patch1:
                endoz1b = EndocyticZone._FromStepsObject(so1, self.mesh)
                self.assertEqual(endoz1b.name, 'endoz1b')
                self.assertEqual(endoz1b.tris, endoz1.tris)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(vesicleGeom))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

