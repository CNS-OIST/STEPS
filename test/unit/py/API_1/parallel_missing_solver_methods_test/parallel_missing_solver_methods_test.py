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

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.utilities.geom_decompose as gd
import steps.mpi
import steps.mpi.solver as ssolver
import steps.rng as srng
from steps.utilities import meshio

class ParallelMissingSolverMethodsTestCase(unittest.TestCase):
    """ 
    Test whether the following methods that were missing from the solver work correctly:
        - setROIAmount
        - setTriCapac
    """
    def setUp(self):
        mdl = smodel.Model()

        S1 = smodel.Spec('S1', mdl)

        vsys = smodel.Volsys('vsys', mdl)
        ssys = smodel.Surfsys('ssys', mdl)

        smodel.Reac('R01', vsys, lhs=[S1], rhs=[S1], kcst=1)
        smodel.SReac('SR01', ssys, slhs=[S1], srhs=[S1], kcst=1)

        vrange = [-200.0e-3, 50e-3, 1e-3]
        vrate = lambda v: 2.0
        Chan1 = smodel.Chan('Chan1', mdl)
        chanop = smodel.ChanState('chanop', mdl, Chan1)
        chancl = smodel.ChanState('chancl', mdl, Chan1)
        smodel.VDepSReac('VDSR01', ssys, slhs=[chancl], srhs=[chanop], k=vrate, vrange=vrange)
        smodel.VDepSReac('VDSR02', ssys, srhs=[chancl], slhs=[chanop], k=vrate, vrange=vrange)

        Chan1_Ohm_I = smodel.OhmicCurr('Chan1_Ohm_I', ssys, chanstate=chanop, g=20e-12, erev=-77e-3)


        if __name__ == "__main__":
            self.mesh = meshio.importAbaqus('meshes/test.inp', 1e-7)[0]
        else:
            self.mesh = meshio.importAbaqus('missing_solver_methods_test/meshes/test.inp', 1e-7)[0]

        comp1 = sgeom.TmComp('comp1', self.mesh, range(self.mesh.countTets()))
        comp1.addVolsys('vsys')

        patch1 = sgeom.TmPatch('patch1', self.mesh, self.mesh.getSurfTris(), comp1)
        patch1.addSurfsys('ssys')

        self.c1ROIInds = range(10)
        self.p1ROIInds = range(5)
        self.mesh.addROI('comp1ROI', sgeom.ELEM_TET, self.c1ROIInds)
        self.mesh.addROI('patch1ROI', sgeom.ELEM_TRI, self.p1ROIInds)

        membrane = sgeom.Memb('membrane', self.mesh, [patch1], opt_method = 1)        

        rng = srng.create('mt19937',512)
        rng.initialize(1234)

        tet_hosts = gd.linearPartition(self.mesh, [1, 1, steps.mpi.nhosts])
        tri_hosts = gd.partitionTris(self.mesh, tet_hosts, self.mesh.getSurfTris())

        self.sim = ssolver.TetOpSplit(mdl, self.mesh, rng, ssolver.EF_DEFAULT, tet_hosts, tri_hosts)
        self.sim.setEfieldDT(1e-4)

        self.sim.reset()
        
    def testSetROIAmount(self):
        amount = 1e-22
        self.assertEqual(self.sim.getROIAmount('comp1ROI', 'S1'), 0)
        self.sim.setROIAmount('comp1ROI', 'S1', amount)
        self.assertAlmostEqual(self.sim.getROIAmount('comp1ROI', 'S1'), amount)

        self.assertEqual(self.sim.getROIAmount('patch1ROI', 'S1'), 0)
        self.sim.setROIAmount('patch1ROI', 'S1', amount)
        self.assertAlmostEqual(self.sim.getROIAmount('patch1ROI', 'S1'), amount)

    def testsetTriCapac(self):
        capac = 1.5e-9
        tri = self.mesh.getSurfTris()[0]
        self.sim.setTriCapac(tri, capac)

    def testgetIClamps(self):
        curr = 1e-12

        tri = self.mesh.getSurfTris()[0]
        self.sim.setTriIClamp(tri, curr)
        self.assertEqual(self.sim.getTriIClamp(tri), curr)

        vert = self.mesh.getTri(tri)[0]
        self.sim.setVertIClamp(vert, curr)
        self.assertEqual(self.sim.getVertIClamp(vert), curr)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ParallelMissingSolverMethodsTestCase, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
