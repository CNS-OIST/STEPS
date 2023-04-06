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

""" Unit tests for steps.geom.Tetmesh class."""

import unittest
from steps.geom import *
import numpy as np

class TetmeshCreationTestCase(unittest.TestCase):
    """ Tests steps.geom.Tetmesh construction. """
    def testCreateTetmesh(self):
        vs=[0.,0.,0.,5.,5.,0.,2.,3.,0.,0.,0.,1.,5.,5.,1.,2.,3.,1.]
        tets=[0,1,2,3,1,4,2,3,2,4,5,3]
        tris=[]

        mesh=Tetmesh(vs,tets,tris)
        self.assertEqual(4*mesh.ntets,len(tets))
        self.assertEqual(3*mesh.nverts,len(vs))


class TetmeshComponentTestCase(unittest.TestCase):
    """ Tests compartment and patch creation. """

    def setUp(self):
        vs=[0.,0.,0.,5.,5.,0.,2.,3.,0.,0.,0.,1.]
        tets=[0,1,2,3]

        self.mesh=Tetmesh(vs,tets,[])

    def testCreateTmComp(self):
        comp=TmComp('comp',self.mesh,range(self.mesh.ntets))
        self.assertEqual(comp.countTets(),self.mesh.ntets)

    def testCreateEmptyTmPatch(self):
        comp=TmComp('comp',self.mesh,range(self.mesh.ntets))
        empty_patch=TmPatch('patch',self.mesh,[],comp)
        self.assertEqual(len(empty_patch.tris),0)

    def testPatchIsTriInside(self):
        comp=TmComp('comp',self.mesh,range(self.mesh.ntets))
        patch12=TmPatch('patch',self.mesh,[1,2],comp)

        self.assertEqual(self.mesh.countTris(),4)
        self.assertEqual([x for x in patch12.tris],[1,2])

        test_tris=[4,1,3,2,1,2,3,4]
        tris_inpatch=[x for x in patch12.isTriInside(test_tris)]
        check_inpatch=[tri in patch12.tris for tri in test_tris]
        self.assertEqual(tris_inpatch,check_inpatch)

    def testIntersection(self):
        vs = [0.,0.,0.,5.,5.,0.,5.,0.,0.,0.,0.,5.]
        tets=[0,1,2,3]
        mesh=Tetmesh(vs, tets, [])

        pts = np.array([[.0, .0, 0.],[.1, .1, 5.1]])
        isecs = mesh.intersect(pts)
        self.assertEqual(len(isecs), 1)
        self.assertEqual(len(isecs[0]), 1)
        self.assertEqual(len(isecs[0][0]), 2)  # tuple
        self.assertEqual(isecs[0][0][0], 0)
        self.assertLess(isecs[0][0][1], 1)
        self.assertGreater(isecs[0][0][1], 0.9)


class TetmeshNPTestCase(unittest.TestCase):
    """ Test numpy-wrapper access with Batch methods. """

    def setUp(self):
        self.vs=[0.,0.,0.,5.,5.,0.,2.,3.,0.,0.,0.,1.,5.,5.,1.,2.,3.,1.]
        self.tets=[0,1,2,3,1,4,2,3,2,4,5,3]
        tris=[]

        self.mesh=Tetmesh(self.vs,self.tets,[])

    def testGetBatchTetsNP(self):
        tets_subset=self.mesh.getBatchTets([0,2])
        self.assertEqual([x for x in tets_subset[0:4]],self.tets[0:4])
        self.assertEqual([x for x in tets_subset[4:8]],self.tets[8:12])

        np_test_subset=np.zeros(8,dtype=np.uint32)
        self.mesh.getBatchTetsNP(np.array([0,2],dtype=np.uint32),np_test_subset)
        np_check=np.array([self.tets[i] for i in [0,1,2,3,8,9,10,11]])
        self.assertEqual([x for x in np_test_subset],[x for x in np_check])



def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TetmeshCreationTestCase, "test"))
    all_tests.append(unittest.makeSuite(TetmeshComponentTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
