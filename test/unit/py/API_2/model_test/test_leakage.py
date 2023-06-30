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
import numpy as np
import os
import math

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *


class Leakage(unittest.TestCase):
    """ Unit tests (based on Rallpack1) for the leakage implementation."""

    # The following values of the analytical solution correspond to the first int(SIM_END/5e-5) values,
    # sampled every int(SAVE_DT/5e-5) steps.  (See analytical solution of the Rallpack1 test)
    # The value 5e-5 is the time step in the analytical solution.
      
    analytic_V_z_min = np.array([-0.065, -0.04247172, -0.03340199, -0.02661633, -0.02103754, -0.0162429, 
                        -0.01201058, -0.0082038, -0.00472958, -0.00152054, 0.00147327, 0.00428992,
                        0.00695834, 0.00950107, 0.01193533, 0.01427481, 0.01653001, 0.01870933,
                        0.02081953, 0.02286587, 0.02485273, 0.02678378, 0.02866195, 0.03048982,
                        0.0322695, 0.03400298])

    analytic_V_z_max = np.array([-0.065, -0.06499991, -0.06496717, -0.06471287, -0.064074, -0.06303987, 
                        -0.06167287, -0.0600496, -0.05823927, -0.05629825, -0.0542707, -0.05219044, 
                        -0.05008313, -0.04796806, -0.0458596, -0.04376848, -0.04170252, -0.03966738, 
                        -0.03766711, -0.03570452, -0.03378142, -0.03189905, -0.03005792, -0.02825834, 
                        -0.02650017, -0.02478309])

    @classmethod
    def setUpClass(cls):
        """Runs the simulation for STEPS 3 and 4. """
        SIM_END = 0.025                # Sim end time (seconds)
        Iinj = 0.1e-9                  # The current injection in amps
        EF_DT = 5e-4
        SAVE_DT = 1e-3
        SEED = 1
        L_G = 0.25                     # Leak conductance, Siemens/m^2
        LEAK_REV = -65.0e-3            # Leak reveral potential, V
        SURFAREA_CYL = 1.0 * math.pi * 1000 * 1e-12
        RA = 1.0                       # Ohm.m
        POT_POS = [0.0, 1.0e-03]

        mesh_path = os.path.join(
                os.path.abspath(__file__ + 5 * '/..'),
                "mesh",
                "axon_cube_L1000um_D866nm_1135tets.msh")

        for USE_STEPS_4 in range(2):
            mesh = DistMesh(mesh_path) if USE_STEPS_4 else TetMesh.LoadGmsh(mesh_path)

            with mesh:
                if USE_STEPS_4:
                    __MESH__ = Compartment.Create()
                    memb = Patch.Create(__MESH__, None)
                    z_min = Patch.Create(__MESH__, None)
                    z_max = Patch.Create(__MESH__, None)
                else:
                    __MESH__ = Compartment.Create(mesh.tets)
                    memb = Patch.Create(mesh.triGroups[(0, "memb")], __MESH__, None)
                    z_min = Patch.Create(mesh.triGroups[(0, "z_min")], __MESH__, None)
                    z_max = Patch.Create(mesh.triGroups[(0, "z_max")], __MESH__, None)


                surfarea_mesh = memb.Area
                corr_fac_area = surfarea_mesh / SURFAREA_CYL

                vol_cyl = math.pi * 0.5 * 0.5 * 1000 * 1e-18
                vol_mesh = __MESH__.Vol
                corr_fac_vol = vol_mesh / vol_cyl

                if USE_STEPS_4:
                    membrane = Membrane.Create([memb], capacitance=0.01 / corr_fac_area)
                    __MESH__.Conductivity = 1 / (RA * corr_fac_vol)
                else:
                    membrane = Membrane.Create([memb])

                POT_TET = TetList(mesh.tets[0, 0, z] for z in POT_POS) # The tetrahedrons from which to record potential
        
            # Model # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
            mdl = Model()
            r = ReactionManager()
            with mdl:
                SA = Species.Create() # STEPS 3 requires at least one species.

            # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
            rng = RNG("r123", 512, SEED)            # Create the solver objects

            if USE_STEPS_4:
                sim = Simulation(
                    "DistTetOpSplit",
                    mdl,
                    mesh,
                    rng,
                    searchMethod=NextEventSearchMethod.GIBSON_BRUCK)
            else:
                part = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)
                sim = Simulation("TetOpSplit", mdl, mesh, rng, MPI.EF_DV_PETSC, part)


            #################### Add leakage ##################
            G_LEAK_SC = corr_fac_area / L_G 
            sim.solver.setMembRes("membrane", G_LEAK_SC, LEAK_REV)

            # Data saving
            rs = ResultSelector(sim)
            Vrs = rs.TETS(POT_TET).V
            sim.toSave(Vrs, dt=SAVE_DT)

            # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
            sim.newRun()
            sim.membrane.Potential = -65e-3

            minzverts = z_min.tris.verts
            sim.VERTS(minzverts).IClamp = Iinj / len(minzverts)

            if not USE_STEPS_4:
                sim.membrane.Capac = 0.01 / corr_fac_area
                sim.membrane.VolRes = RA * corr_fac_vol

            sim.EfieldDT = EF_DT
            sim.run(SIM_END)

            cls.sim = sim
            if not USE_STEPS_4:
                cls.Vzmin3 = Vrs.data[0, :, 0] 
                cls.Vzmax3 = Vrs.data[0, :, 1]
            else:
                cls.Vzmin4 = Vrs.data[0, :, 0] 
                cls.Vzmax4 = Vrs.data[0, :, 1]

    def testMagnitude(self):
        """Test getter and setter of membrane resistance in STEPS4"""
        res = 0.4
        rev_pot = 0.4
        self.sim.solver.setMembRes("membrane", res, rev_pot)
        self.assertEqual( (res, rev_pot), self.sim.solver.getMembRes("membrane"), 
                        "Error in set or get method")

    def testSolution_not_constant(self):
        self.assertTrue( np.any( np.abs(self.Vzmax3-self.Vzmax3[0]) > 1e-15  ), 
                        "The solution is constant")

    def testMse_3vs4(self):
        """"MSE between the values obtained from steps 3 and steps 4"""
        threshold = 1e-12
        MSE_max = np.mean( (self.Vzmax3 - self.Vzmax4)**2 )
        self.assertLess(MSE_max , threshold, "Error MSE(Vzmax3-Vzmax4) bigger than threshold")
        MSE_min = np.mean( (self.Vzmin3 - self.Vzmin4)**2 )
        self.assertLess(MSE_min , threshold, "Error MSE(Vzmin3-Vzmin4) bigger than threshold")
        self.assertNotEqual(MSE_min, 0.0, "MSE equal zero")
        self.assertNotEqual(MSE_max, 0.0, "MSE equal zero")

    def testMse_analytical(self):
        """"MSE between the values obtained from steps 3 steps 4 and the analytical solution"""
        threshold = 1e-2

        self.assertEqual( len(self.analytic_V_z_max), len(self.Vzmax3), 
                        "analytic and numeric solution have different length")

        MSE_max = np.mean( (self.Vzmax3 - self.analytic_V_z_max)**2 )
        self.assertLess(MSE_max , threshold, "Error MSE(Vzmax3-analytic) bigger than threshold")
        MSE_min = np.mean( (self.Vzmin3 - self.analytic_V_z_min)**2 )
        self.assertLess(MSE_min , threshold, "Error MSE(Vzmin3-analytic) bigger than threshold")

        MSE_max = np.mean( (self.Vzmax4 - self.analytic_V_z_max)**2 )
        self.assertLess(MSE_max , threshold, "Error MSE(Vzmax4-analytic) bigger than threshold")
        MSE_min = np.mean( (self.Vzmin4 - self.analytic_V_z_min)**2 )
        self.assertLess(MSE_min , threshold, "Error MSE(Vzmin4-analytic) bigger than threshold")


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(Leakage, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
