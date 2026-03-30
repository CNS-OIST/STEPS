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

import argparse
import importlib
import os
import sys
import tempfile
import unittest

import mpi4py
import numpy as np
import steps.interface
from steps.saving import *
from steps.sim import *

from steps.geom import *
from steps.model import *
from steps.rng import *

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.join(LOCAL_DIR, "..", "..", "..")
MESH_DIR = os.path.join(TEST_DIR, "mesh")

NA = 6.0221415e23


def getUniqueTempPrefix(*args, **kwargs):
    _, pathPrefix = tempfile.mkstemp(*args, **kwargs)
    os.remove(pathPrefix)
    if MPI._nhosts > 1:
        # Need to synchronize across ranks
        pathPrefix = mpi4py.MPI.COMM_WORLD.bcast(pathPrefix, root=0)
    return pathPrefix


def volDiffFunc(NINJECT, DCST, r, t):
    Crt = NINJECT * np.exp(-(r**2) / (4 * DCST * t)) / (8 * (np.pi * DCST * t) ** 1.5)
    # Convert from number of molecules per m^3 to M
    return Crt * (1e-3 / NA)


def firstOrderIrrevFunc(C0, KCST, t):
    return C0 * np.exp(-KCST * t)


def diffTubeInfluxFunc(DCST, KCST, N, Area, L, x, t):
    J = KCST * N / Area * 1e-3
    nmax = 100
    sumvalue = sum(
        (-1) ** n
        / (n**2)
        * np.exp(-DCST * t * (n * np.pi / L) ** 2)
        * np.cos(n * np.pi * x / L)
        for n in range(1, nmax)
    )
    return (
        J
        * L
        / DCST
        * (
            DCST * t / (L**2)
            + (3 * x**2 - L**2) / (6 * L**2)
            - 2 / (np.pi**2) * sumvalue
        )
    )


def diffTubeDegradationFunc(C0, DCST, L, d, time):
    d = abs(d)
    nmax = 100
    sumvalue = sum(
        1
        / (2 * n + 1)
        * np.exp(-(DCST * (2 * n + 1) ** 2 * np.pi**2 * time) / (4 * L**2))
        * np.sin((2 * n + 1) * np.pi * d / (2 * L))
        for n in range(nmax)
    )
    return 4 * C0 / np.pi * sumvalue


def get_average_concentrations(concs, tind, nbins):
    dists = concs.metaData["distance"]
    vols = np.array(concs.metaData["volume"])
    bins = np.histogram_bin_edges(dists, bins=nbins)
    grps = np.digitize(dists, bins)
    counts = np.bincount(grps)[1:-1]

    avg_pos = np.bincount(grps, weights=dists)[1:-1] / counts
    total_vol = np.bincount(grps, weights=vols)[1:-1]
    avg_conc = (
        np.bincount(grps, weights=np.mean(concs.data[:, tind, :] * vols, axis=0))[1:-1]
        / total_vol
    )

    return avg_pos, avg_conc


class DistTetopsplitRLeaping(unittest.TestCase):
    """Tests the accuracy of the R-leaping and diff-leaping methods for DistTetOpSplit"""

    def setUp(self):
        sys.path.append(LOCAL_DIR)
        self.createdFiles = set()

    def tearDown(self):
        if MPI.rank == 0:
            for f in self.createdFiles:
                if os.path.exists(f):
                    os.remove(f)

    def get_common_args(self):
        args = argparse.Namespace()
        args.ENDT = 0.1
        args.DT = 0.01
        args.NITER = 100
        args.seed = 4
        args.noscaletime = False
        args.meshPath = os.path.join(MESH_DIR, "box.msh")
        args.meshScale = 1e-6
        args.indepKProcs = False
        args.diffTolerance = 0.1
        args.skellamThresh = 1.0
        args.maxskips = 5
        args.diffmindtfactor = 1.0
        args.diffCrankNicolson = 2.0
        args.ssaThreshold = 50
        args.reacTolerance = 0.05
        args.reacTheta = 0.1
        return args

    def run_simulations(self, model_name, args):
        """Run simulations of the given model"""
        model = importlib.import_module(f"models.{model_name}")

        ########################################################################

        params = dict(DT=args.DT, ENDT=args.ENDT, timefact=1)

        model.set_parameters(args, params)

        mdl = Model()
        model.setup_model(mdl, args, params)

        ########################################################################

        params["meshPath"] = args.meshPath
        params["meshScale"] = args.meshScale
        mesh = DistMesh(params["meshPath"], params["meshScale"])

        model.setup_mesh(mesh, mdl, args, params)

        ########################################################################

        rng = RNG("mt19937", 512, args.seed)

        sim = Simulation(
            "DistTetOpSplit",
            mdl,
            mesh,
            rng,
            SSAMethod=SSAMethod.RLEAPING,
            searchMethod=NextEventSearchMethod.RLEAPING,
            diffMethod=DiffusionMethod.TAU_LEAPING_DT,
            # diffMethod=DiffusionMethod.CONSTANT_DT,
            check=False,
            indepKProcs=args.indepKProcs,
        )
        sim.setDiffusionTolerance(args.diffTolerance)
        sim.setDiffusionMaxDtSkips(args.maxskips)
        sim.setDiffusionMinDtFactor(args.diffmindtfactor)
        sim.setDiffusionNormalApproximationThreshold(args.skellamThresh)
        sim.setDiffusionCrankNicolsonThreshold(args.diffCrankNicolson)
        sim.setReactionSSAThreshold(args.ssaThreshold)
        sim.setReactionTolerance(args.reacTolerance)
        sim.setReactionTheta(args.reacTheta)

        ########################################################################

        model.setup_sim(sim, args, params)

        ########################################################################

        pathPrefix = getUniqueTempPrefix(prefix=f"{self.__class__.__name__}Validation")
        self.createdFiles.add(f"{pathPrefix}.h5")

        with HDF5Handler(pathPrefix, mode="w") as hdf:
            group = sim.toDB(hdf, "validations", **vars(args))

            model.setup_group(group, sim, args, params)

            for i in range(args.NITER):
                sim.newRun()
                model.initialize_sim(sim, args, params)
                sim.run(params["ENDT"])

        # Make sure that all ranks finished writing their files
        mpi4py.MPI.COMM_WORLD.barrier()

        return pathPrefix

    def run_and_compute_relative_error(
        self, model, args, tind=-1, nbins=30, do_plot=False
    ):
        pathPrefix = self.run_simulations(model, args)

        if MPI.rank == 0:
            # Only rank 0 reads data
            with HDF5Handler(pathPrefix) as hdf:
                group = hdf.get()
                concs = group.results["Concentrations"]

                time = concs.time[0, -1]
                INJ = group.parameters.get("INJ", 0)
                DCST = group.parameters.get("DCST", 1)
                KCST = group.parameters.get("KCST", 1)
                volumes = concs.metaData["volume"]
                VOL = np.sum(volumes)
                match model:
                    case "DiffPoint":
                        xvals, experimental = get_average_concentrations(
                            concs, tind, nbins
                        )
                        analytical = np.array(
                            [volDiffFunc(INJ, DCST, d, time) for d in xvals]
                        )
                    case "DiffUniform":
                        ntets = len(concs.data[0, tind, :])
                        nruns = len(concs.data[:, tind, 0])
                        mu = np.mean(concs.data[:, tind, :], axis=0)
                        sigma = np.sqrt(nruns / (nruns - 1)) * np.std(
                            concs.data[:, tind, :], axis=0
                        )
                        experimental = np.mean(sigma[mu != 0] / mu[mu != 0])
                        analytical = np.mean(
                            [
                                np.sqrt(INJ * (v / VOL)) / (INJ * (v / VOL))
                                for v in volumes
                            ]
                        )
                        if do_plot:
                            from matplotlib import pyplot as plt

                            plt.bar([1], [experimental], label="simulation")
                            plt.axhline(analytical, color="r", label="analytical")
                            plt.legend()
                            plt.show()
                        return mpi4py.MPI.COMM_WORLD.bcast(
                            abs(experimental - analytical) / analytical, root=0
                        )
                    case "ReacFirstOrderIrrev":
                        xvals = concs.time[0]
                        experimental = np.mean(
                            np.sum(concs.data[:, :, :] * volumes, axis=2)
                            / np.sum(volumes),
                            axis=0,
                        )
                        analytical = np.array(
                            [
                                firstOrderIrrevFunc(INJ / NA / VOL * 1e-3, KCST, t)
                                for t in xvals
                            ]
                        )
                    case "DiffTubeInflux":
                        Area = group.staticData["BorderArea"]
                        ntets = group.staticData["BorderNbTets"]
                        L = group.staticData["Length"] / 2
                        xvals, experimental = get_average_concentrations(
                            concs, tind, nbins
                        )
                        analytical = np.array(
                            [
                                diffTubeInfluxFunc(
                                    DCST, KCST / NA, ntets, Area, L, d, time
                                )
                                for d in xvals
                            ]
                        )
                    case "DiffTubeDegradation":
                        L = group.staticData["Length"] / 2
                        xvals, experimental = get_average_concentrations(
                            concs, tind, nbins
                        )
                        analytical = np.array(
                            [
                                diffTubeDegradationFunc(
                                    2 * INJ / NA / VOL * 1e-3, DCST, L, d, time
                                )
                                for d in xvals
                            ]
                        )

                if do_plot:
                    from matplotlib import pyplot as plt

                    plt.plot(xvals, experimental, "--", label="simulation")
                    plt.plot(xvals, analytical, "r-", label="analytical")
                    plt.legend()
                    plt.show()

                return mpi4py.MPI.COMM_WORLD.bcast(
                    np.sqrt(np.mean((experimental - analytical) ** 2))
                    / np.mean(analytical),
                    root=0,
                )
        else:
            return mpi4py.MPI.COMM_WORLD.bcast(0, root=0)

    def test_diffusion_point(self):
        """Test that point diffusion works as expected"""
        args = self.get_common_args()
        args.DCST = 1e-10
        args.ENDT = 0.05
        args.DT = 0.05
        args.meshScale = 1e-5

        injs = [100, 1e4, 1e6]
        # We need more iterations when there are fewer molecules injected
        iters = [1000, 50, 1]
        for inj, iter in zip(injs, iters):
            args.INJ = inj
            args.NITER = iter
            with self.subTest(INJ=inj, NITER=iter):
                error = self.run_and_compute_relative_error("DiffPoint", args)
                if MPI.rank == 0:
                    print(f"{inj=}, {iter=}, {error=}")
                self.assertLessEqual(error, 0.1)

    def test_diffusion_uniform(self):
        """Inject species in the whole mesh and let them diffuse"""
        args = self.get_common_args()
        args.DCST = 1e-10
        args.ENDT = 0.05
        args.DT = 0.05
        args.meshScale = 1e-5

        injs = [100, 1e4, 1e6, 1e9]
        # We need more iterations when there are fewer molecules injected
        iters = [1000, 10, 10, 10]
        for inj, iter in zip(injs, iters):
            args.INJ = inj
            args.NITER = iter
            with self.subTest(INJ=inj, NITER=iter):
                error = self.run_and_compute_relative_error("DiffUniform", args)
                if MPI.rank == 0:
                    print(f"{inj=}, {iter=}, {error=}")
                self.assertLessEqual(error, 0.1)

    def test_diffusion_tube_degradation(self):
        """Two separated reactants diffuse towards each other in the center of a cylinder and annihilate."""
        args = self.get_common_args()
        args.DCST = 1e-10
        args.ENDT = 0.05
        args.DT = 0.05
        args.NITER = 10
        args.INJ = 1e6
        args.meshPath = os.path.join(MESH_DIR, "10um_tube_4714tets.msh")
        args.meshScale = 1e-6

        error = self.run_and_compute_relative_error("DiffTubeDegradation", args)
        if MPI.rank == 0:
            print(f"{error=}")
        self.assertLessEqual(error, 0.02)

    def test_tube_influx(self):
        """1D diffusion in a finite tube with a constant and equal influx of the same species
        of molecule at both ends"""
        args = self.get_common_args()
        args.KCST = 1000
        args.DCST = 1e-10
        args.ENDT = 0.5
        args.DT = 0.5
        args.NITER = 10
        args.INJ = 1e6
        args.meshPath = os.path.join(MESH_DIR, "10um_tube_4714tets.msh")
        args.meshScale = 1e-6

        error = self.run_and_compute_relative_error("DiffTubeInflux", args)
        if MPI.rank == 0:
            print(f"{error=}")
        self.assertLessEqual(error, 0.05)

    def test_reac_first_order_irrev(self):
        """Inject species in the whole mesh and let them degrade through a first order reaction.
        Diffusion is present but does not play much of a role."""
        args = self.get_common_args()
        args.KCST = 100
        args.DCST = 1e-11
        args.ENDT = 0.1
        args.DT = 0.01
        args.NITER = 10
        args.INJ = 1e6
        args.meshScale = 1e-5

        error = self.run_and_compute_relative_error("ReacFirstOrderIrrev", args)
        if MPI.rank == 0:
            print(f"{error=}")
        self.assertLessEqual(error, 0.01)


def suite():
    all_tests = []
    all_tests.append(
        unittest.TestLoader().loadTestsFromTestCase(DistTetopsplitRLeaping)
    )
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
