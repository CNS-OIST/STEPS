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

import os.path as path
import pathlib
import unittest

import numpy as np
import steps.interface
from steps.sim import *

from . import rallpack3


class TestRallpack3Base(unittest.TestCase):

    V0_GHK_REF = "rallpack3_ghk_0_0.001dt_1000seg"
    V1_GHK_REF = "rallpack3_ghk_1000_0.001dt_1000seg"

    # defaults
    default_params = {
        "plot": False,
        "meshdir": "validation_efield/meshes",
        "mesh": "axon_ecs_1000um_1um_1314tets.msh",
        "datadir": "validation_efield/data/rallpack3_benchmark/",
        "v0data": "rallpack3_0_0.001dt_1000seg",
        "v1data": "rallpack3_1000_0.001dt_1000seg",
        "seed": 10,
        "saving_mode": "none",
        "checks": ["start", "end"],
        "max_rms_err": 2,
    }

    def _test_rallpack3(self, **params):
        params = {**self.default_params, **params}

        meshfile = path.join(params["meshdir"], params["mesh"])

        time, Vzmin, Vzmax, V0ref, V1ref = rallpack3.run_comparison(meshfile, params)

        rms_err_0um = np.linalg.norm(Vzmin[0, :, 0] - V0ref) / np.sqrt(V0ref.shape[0])
        rms_err_1000um = np.linalg.norm(Vzmax[0, :, 0] - V1ref) / np.sqrt(
            V1ref.shape[0]
        )

        print(f"rms error at 0um = {rms_err_0um}")
        print(f"rms error at 1000um = {rms_err_1000um}")

        if params["plot"] and MPI.rank == 0:
            import matplotlib.pyplot as plt

            plt.subplot(211)
            plt.plot(time[0, :], V0ref, "k-", label="Correct, 0um", linewidth=3)
            plt.plot(
                time[0, :],
                Vzmin[0, :, 0],
                "r--",
                label="STEPS, 0um",
                linewidth=3,
            )
            plt.legend(loc="best")
            plt.ylabel("Potential (mV)")
            plt.subplot(212)
            plt.plot(time[0, :], V1ref, "k-", label="Correct, 1000um", linewidth=3)
            plt.plot(
                time[0, :],
                Vzmax[0, :, 0],
                "r--",
                label="STEPS, 1000um",
                linewidth=3,
            )
            plt.legend(loc="best")
            plt.ylabel("Potential (mV)")
            plt.xlabel("Time (ms)")
            plt.show(block=True)

        if "start" in params["checks"]:
            self.assertLess(rms_err_0um, params["max_rms_err"])
        if "end" in params["checks"]:
            self.assertLess(rms_err_1000um, params["max_rms_err"])

        return time, Vzmin, Vzmax


class TestRallpack3(TestRallpack3Base):

    def test_rallpack3_tetexact(self):
        self._test_rallpack3(solver="Tetexact", currType="Ohmic", plot=False)

    # This test can generate the reference data for GHK currents (with write_data=True)
    def test_rallpack3_tetODE(self, write_data=False):
        time, Vzmin, Vzmax = self._test_rallpack3(
            solver="TetODE",
            currType="GHK",
            SAVE_DT=5e-6,
            ion_diff=True,
            plot=False,
            v0data=self.V0_GHK_REF,
            v1data=self.V1_GHK_REF,
            max_rms_err=1,
        )

        if write_data:
            datadir = pathlib.Path(self.default_params["datadir"])
            paths = [datadir.joinpath(self.V0_GHK_REF), datadir.joinpath(self.V1_GHK_REF)]
            data = [Vzmin, Vzmax]
            for i, path, vals in zip(range(2), paths, data):
                with path.open("w") as f:
                    f.write(
                        "\n".join(
                            [f"label:v({i})", f"{vals.shape[1]}"]
                            + [f"{t} {v}" for t, v in zip(time[0, :], vals[0, :, 0])]
                        )
                    )

    # Disabled because too slow
    # def test_rallpack3_ghk_tetexact(self):
    #     self._test_rallpack3(solver="Tetexact", currType="GHK", plot=True, max_rms_err=15)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestRallpack3))
    return unittest.TestSuite(all_tests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
