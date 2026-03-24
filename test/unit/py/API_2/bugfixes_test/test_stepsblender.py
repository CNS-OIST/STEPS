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

"""Bugfixes test for stepsblender package"""

import tempfile
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.rng import *
from steps.saving import *


import os

FILEDIR = os.path.dirname(os.path.abspath(__file__))
MESHDIR = os.path.join(FILEDIR, '..', '..', '..', '..', 'mesh')

DT = 0.01
END_T = 0.1

class StepsBlenderBugFixes(unittest.TestCase):
    """Tests related to bugfixes in the stepsblender package"""
    def setUp(self):
        self.createdFiles = set()

    def tearDown(self):
        for path in self.createdFiles:
            if os.path.isfile(path):
                os.remove(path)

    def testUnattributedTetsBugFix(self):
        """Test that stepsblender does not crash when loading data that contains tetrahedrons not associated with any compartment"""
        try:
            from stepsblender.dataloader import _HDF5BlenderDataLoader
        except ImportError:
            return

        mdl = Model()
        with mdl:
            vsys = VolumeSystem.Create()
            SA = Species.Create()
            with vsys:
                Diffusion(SA, 0)
            
        mesh = TetMesh.LoadGmsh(os.path.join(MESHDIR, 'cube.msh'), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets[:-2])

        rng = RNG()
        sim = Simulation('Tetexact', mdl, mesh, rng)

        rs = ResultSelector(sim)

        sim.toSave(rs.TETS().SA.Count, dt=DT)

        _, hdfPath = tempfile.mkstemp(prefix=f'{self.__class__.__name__}teststepsblender', suffix='')
        self.createdFiles.add(hdfPath)
        with XDMFHandler(hdfPath) as hdf:
            sim.toDB(hdf)
            sim.newRun()
            sim.run(END_T)
            self.createdFiles |= set(hdf._getFilePaths())

        # Load data with stepsblender
        with HDF5Handler(hdfPath) as hdf:
            group = hdf.get()
            meshGroup = group._group[XDMFHandler._MESH_GROUP_NAME]
            loader = _HDF5BlenderDataLoader(hdfPath, None, None)
            loader._getMesh(hdf, group, meshGroup, False)

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(StepsBlenderBugFixes))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

