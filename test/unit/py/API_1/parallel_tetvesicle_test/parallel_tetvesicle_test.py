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



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.mpi
import steps.mpi.solver as solv
from steps.utilities import meshio
import steps.utilities.geom_decompose as gd

import numpy as np

class ParallelTetVesicleTestCase(unittest.TestCase):
    """ Test cases for parallel TetVesicle solver. """

    def setUp(self):
        if __name__ == "__main__":
            self.mesh = meshio.loadMesh('meshes/axon_cube_L1000um_D866nm_1978tets')[0]
        else:
            self.mesh = meshio.loadMesh('parallel_tetvesicle_test/meshes/axon_cube_L1000um_D866nm_1978tets')[0]

        self.param = {
        ## Rallpack1:
            'R_A'  :    1.0,          # axial resistivity Ω·m
            'R_M'  :    4.0,          # membrane resistivity Ω·m²
            'C_M'  :    0.01,         # membrane capacity F/m²
            'E_M'  :   -0.065,        # p.d. across membrane V
            'Iinj' :    0.1e-9,       # injection current A
            'diameter': 1.0e-6,       # cylinder diameter m
            'length':   1.0e-3,       # cylinder length m
        # STEPS
            'sim_end':  0.25,         # simulation stop time s
            'EF_dt':    5.0e-5,       # E-field evaluation time step s
        }

        cyto = sgeom.TmComp('cyto', self.mesh, range(self.mesh.ntets))
        memb_tris = list(self.mesh.getSurfTris())
        memb = sgeom.TmPatch('memb', self.mesh, memb_tris, cyto)
        membrane = sgeom.Memb('membrane', self.mesh, [memb] )

        self.mdl = smodel.Model()
        memb = sgeom.castToTmPatch(self.mesh.getPatch('memb'))

        ssys = smodel.Surfsys('ssys', self.mdl)
        memb.addSurfsys('ssys')

        L = smodel.Chan('L', self.mdl)
        Leak = smodel.ChanState('Leak', self.mdl, L)

        # membrane conductance
        area_cylinder = np.pi * self.param['diameter'] * self.param['length']
        L_G_tot = area_cylinder / self.param['R_M']
        g_leak_sc = L_G_tot / len(memb.tris)
        OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate = Leak, erev = self.param['E_M'], g = g_leak_sc)

        self.rng = srng.create('r123', 512)
        self.rng.initialize(1000)

    def tearDown(self):
        self.mdl = None
        self.mesh = None
        self.rng = None

    def testCreateSolver(self):
        """
        Try to create and initialize TetVesicle solver.
        """
        solver = solv.TetVesicle(self.mdl, self.mesh, self.rng, solv.EF_DV_PETSC)
        if steps.mpi.rank == 0:
            self.assertEqual(solver.getSolverName(), 'TetVesicleVesRaft')
        else:
            self.assertEqual(solver.getSolverName(), 'TetVesicleRDEF')


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(ParallelTetVesicleTestCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
