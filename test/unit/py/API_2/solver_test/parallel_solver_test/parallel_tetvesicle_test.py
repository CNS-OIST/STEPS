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

""" Unit tests for parallel tetvesicle solver."""

import os
import unittest

import numpy as np

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

class ParallelTetVesicleTestCase(unittest.TestCase):
    """Test parallel TetVesicle."""

    def setUp(self):
        dirPath = os.path.dirname(os.path.abspath(__file__))
        meshPath = os.path.join(dirPath, 'meshes/axon_cube_L1000um_D866nm_1978tets')

        self.mesh = TetMesh.Load(meshPath)

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
        # membrane conductance
        area_cylinder = np.pi * self.param['diameter'] * self.param['length']
        L_G_tot = area_cylinder / self.param['R_M']
        g_leak_sc = L_G_tot / len(self.mesh.surface)

        self.mdl = Model()
        with self.mdl:
            ssys = SurfaceSystem.Create()

            Leak = SubUnitState.Create()
            L = Channel.Create([Leak])

            with ssys:
                OC_L = OhmicCurr.Create(L[Leak], g_leak_sc, self.param['E_M'])
        
        with self.mesh:
            cyto = Compartment.Create(self.mesh.tets)
            memb = Patch.Create(self.mesh.surface, cyto, None, ssys)

            membrane = Membrane.Create([memb])

        self.rng = RNG('r123', 512, 1000)

    def testCreateSolver(self):
        """Try to create and initialize TetVesicle solver"""

        sim = Simulation('TetVesicle', self.mdl, self.mesh, self.rng, MPI.EF_DV_PETSC)
        if MPI.rank == 0:
            self.assertEqual(sim.solver.getSolverName(), 'TetVesicleVesRaft')
        else:
            self.assertEqual(sim.solver.getSolverName(), 'TetVesicleRDEF')

def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(ParallelTetVesicleTestCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

