########################################################################

# Stochastic degradation-diffusion process.

# AIMS: to verify if zero order reactions are inserted correctly to the SSA
# after bugfix https://github.com/CNS-OIST/HBP_STEPS/pull/146
# Note: 
# 1. This is not validated against analytical result, but 
# only check if the zero order reaction can be found through SSA search.
# Future improvement is needed.
# 2. sim.reset() should be avoided in this test
  
########################################################################

import unittest

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import math
import time

########################################################################

class TestRDWellMixedWmrk4(unittest.TestCase):
    def get_model(self):
        mdl = Model()
        r = ReactionManager()
        with mdl:
            SA = Species.Create()
            volsys = VolumeSystem.Create()
            with volsys:
                None >r['R1']> SA
                r['R1'].K = 0.001
        return mdl

    def get_geom(self, geom_type):
        if geom_type == "well_mixed":
            geom = Geometry()
            with geom:
                comp = Compartment.Create('volsys', 1e-18)
            return geom
        elif geom_type == "tetmesh":
            mesh_loc = None
            if __name__== "__main__":
                mesh_loc = "meshes/1x1x1.inp"
            else:
                mesh_loc = "validation_rd/meshes/1x1x1.inp"
            mesh = TetMesh.LoadAbaqus(mesh_loc, scale=1e-06)
            with mesh:
                comp = Compartment.Create(mesh.tets, 'volsys')
            return mesh
        else:
            self.assertTrue(False)

    def get_solver(self, solver_type, mdl, geom):
        r = RNG('r123', 1000, 1234)
        if solver_type == "Wmrk4":
            solver_Modif = Simulation('Wmrk4', mdl, geom)
            solver_Modif.setDT(1e-5)
            return solver_Modif
        elif solver_type == "Wmdirect":
            solver_Modif = Simulation('Wmdirect', mdl, geom, r)
            return solver_Modif
        elif solver_type == "Wmrssa":
            solver_Modif = Simulation('Wmrssa', mdl, geom, r)
            return solver_Modif
        elif solver_type == "Tetexact":
            solver_Modif = Simulation('Tetexact', mdl, geom, r)
            return solver_Modif
        else:
            self.assertTrue(False)

    def validate(self, sim):
        sim.newRun()
        sim.run(1)
        count = sim.comp.SA.Count
        self.assertTrue(count > 0)
        if __name__== "__main__":
            print("Count: ", count)

    def test_zeroorder_wmrk4(self):
        "Zero order reaction (Wmrk4)"
        if __name__== "__main__":
            print("Zero order reaction (Wmrk4)")
        mdl = self.get_model()
        geom = self.get_geom("well_mixed")
        sim = self.get_solver("Wmrk4", mdl, geom)
        self.validate(sim)

    def test_zeroorder_wmdirect(self):
        "Zero order reaction (Wmdirect)"
        if __name__== "__main__":
            print("Zero order reaction (Wmdirect)")
        mdl = self.get_model()
        geom = self.get_geom("well_mixed")
        sim = self.get_solver("Wmdirect", mdl, geom)
        self.validate(sim)

    def test_zeroorder_wmrssa(self):
        "Zero order reaction (Wmrssa)"
        if __name__== "__main__":
            print("Zero order reaction (Wmrssa)")
        mdl = self.get_model()
        geom = self.get_geom("well_mixed")
        sim = self.get_solver("Wmrssa", mdl, geom)
        self.validate(sim)

    def test_zeroorder_tetexact(self):
        "Zero order reaction (Tetexact)"
        if __name__== "__main__":
            print("Zero order reaction (Tetexact)")
        mdl = self.get_model()
        geom = self.get_geom("tetmesh")
        sim = self.get_solver("Tetexact", mdl, geom)
        self.validate(sim)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestRDWellMixedWmrk4, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
