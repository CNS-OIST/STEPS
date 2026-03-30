import unittest
 
import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import numpy
import math
import time
import os
import mpi4py.MPI as MPI


class TestTetWeighting(unittest.TestCase):

    def setUp(self):

        # Simple model with a few reactions and diffusions
        TF = 10
        CF = 5e4
        DCST = 1e-13
        A = 0.3
        B = 3
        self.model = Model()
        r = ReactionManager()
        with self.model:
            X, Y = Species.Create()
            vsys = VolumeSystem.Create()
            with vsys:
                None < r[1] > X > r[2] > Y
                r[1].K = A * TF / CF, TF
                r[2].K = B * TF

                2*X + Y > r[3]> 3*X
                r[3].K = TF * (CF ** 2)

                Diffusion(X, DCST)
                Diffusion(Y, DCST)
        
        self.mesh = TetMesh.LoadAbaqus('validation_rd/meshes/brick_40_4_4_1400tets.inp', scale=1e-06)
        with self.mesh:
            comp = Compartment.Create(self.mesh.tets, vsys)

        self.rng = RNG('mt19937', 512, 1234)

    def test_extent_weighting(self):
        SIM_TYPES = ['TetOpSplit', 'TetVesicle']
        ENDT = 1.5

        # Test for each simulation type
        for sim_type in SIM_TYPES:
            print(f'Running extent weighting test for {sim_type}...')

            ### Initial collection run; gather weights ###
            sim = Simulation(sim_type, self.model, self.mesh, self.rng)
            sim.newRun()
            sim.comp.X.Conc = 5e-6
            sim.comp.Y.Conc = 1.6e-4
            for j in range(1000):
                t = ENDT * j / 999
                sim.run(t)
            sim.saveTetWeights(prefix='validation_tet_weighting/')
            MPI.COMM_WORLD.Barrier()  # Ensure file is written before proceeding

            ### Load gathered weights, use weighted mesh, run again ###
            start_host = 0 if sim_type == 'TetOpSplit' else 1
            partition = TetWeightPartition(self.mesh, prefix="validation_tet_weighting/", solver=sim_type)

            # Check that partition can be loaded in a new simulation
            sim2 = Simulation(sim_type, self.model, self.mesh, self.rng, partition)

            ### Assertions for reasonable weighted partition ###
            # Load weight data according to partitioning
            host_range = range(0, MPI.COMM_WORLD.Get_size())
            weights_per_host = {i:0 for i in host_range}
            with open('validation_tet_weighting/weights.log') as f:
                for idx, weight in enumerate(f):
                    host = partition.tetPart[idx]
                    weights_per_host[host] += float(weight)

            # If TetVesicle, host 0 should have zero weight
            if start_host == 1:
                self.assertEqual(weights_per_host[0], 0.0)

            # Check that total weights between partitions are reasonably balanced
            ideal_division = 1 / (MPI.COMM_WORLD.Get_size() - start_host)
            deviation = ideal_division * 0.2
            lwb = ideal_division - deviation
            upb = ideal_division + deviation
            total_weight = sum(weights_per_host.values())
            for host, w in weights_per_host.items():
                if host < start_host:
                    continue
                self.assertTrue(lwb < w / total_weight < upb, f"Host {host} weight fraction {w / total_weight} outside bounds [{lwb}, {upb}]")

            ### Remove weight log file ###
            # Make sure all processes done using file before removing
            MPI.COMM_WORLD.Barrier()
            if MPI.COMM_WORLD.Get_rank() == 0:
                os.remove('validation_tet_weighting/weights.log')
            MPI.COMM_WORLD.Barrier()


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestTetWeighting))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
