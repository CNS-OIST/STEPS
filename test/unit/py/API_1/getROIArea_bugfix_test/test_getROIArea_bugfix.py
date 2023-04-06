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

# Unit test for bugfix https://github.com/CNS-OIST/HBP_STEPS/pull/223

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.model as smodel
import steps.geom as sgeom
from steps.utilities import meshio

class getROIAreaTestCase(unittest.TestCase):
    """ Test if the getROIArea method returns the area without throwing an error. """
    def setUp(self):
        self.model = smodel.Model()

        self.vsys1 = smodel.Volsys('vsys1', self.model)
        self.ssys1 = smodel.Surfsys('ssys1', self.model)
    
        if __name__ == "__main__":
            self.mesh = meshio.loadMesh('meshes/cyl_len10_diam1')[0]
        else:
            self.mesh = meshio.loadMesh('getROIArea_bugfix_test/meshes/cyl_len10_diam1')[0]

        ntets = self.mesh.countTets()
        comp1Tets, comp2Tets = [], []
        comp1Tris, comp2Tris = set(), set()
        for i in range(ntets):
            if self.mesh.getTetBarycenter(i)[0] > 0:
                comp1Tets.append(i)
                comp1Tris |= set(self.mesh.getTetTriNeighb(i))
            else:
                comp2Tets.append(i)
                comp2Tris |= set(self.mesh.getTetTriNeighb(i))
        patch1Tris = list(comp1Tris & comp2Tris)

        self.comp1 = sgeom.TmComp('comp1', self.mesh, comp1Tets)
        self.comp2 = sgeom.TmComp('comp2', self.mesh, comp2Tets)
        self.comp1.addVolsys('vsys1')
        self.comp2.addVolsys('vsys1')

        self.patch1 = sgeom.TmPatch('patch1', self.mesh, patch1Tris, self.comp1, self.comp2)
        self.patch1.addSurfsys('ssys1')

        self.ROI1 = self.mesh.addROI('ROI1', sgeom.ELEM_TRI, patch1Tris)
    
    def tearDown(self):
        self.model = None
        self.mesh = None
    
    def testgetROIArea(self):
        area = self.mesh.getROIArea('ROI1')

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(getROIAreaTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
