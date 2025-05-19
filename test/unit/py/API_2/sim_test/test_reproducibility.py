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

""" Unit tests for tetexact reproducibility."""

import importlib
import numpy as np
import os
import tempfile
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

from . import base_model


class ReproducibilityTestCase:
    """Test reproducibility of STEPS solvers"""

    def setUp(self):
        super().setUp()

        self.nbRuns = 5
        self.createdFiles = set()

    def tearDown(self):
        super().tearDown()
        if MPI._shouldWrite:
            for path in self.createdFiles:
                if os.path.isfile(path):
                    os.remove(path)

    def _getResultsSelector(self, rs):
        return rs.ALL(Compartment, Patch).ALL(Species).Count

    def _testReproducible(self, plot=False):
        _, hdfprefix = tempfile.mkstemp(prefix=f'{self.__class__.__name__}data')
        hdfpath = hdfprefix + '.h5'
        self.createdFiles.add(hdfpath)
        if MPI._shouldWrite and os.path.exists(hdfpath):
            os.remove(hdfpath)

        for i in range(self.nbRuns):
            mdl = self.get_API2_Mdl()
            geom = self.get_API2_Geom(mdl)
            sim = self._get_API2_Sim(mdl, geom)

            rs = ResultSelector(sim)
            res = self._getResultsSelector(rs)

            sim.toSave(res, dt=self.deltaT)

            with HDF5Handler(hdfprefix) as hdf:
                sim.toDB(hdf, 'test')
                sim.newRun()
                self.init_API2_sim(sim)

                sim.run(self.endTime / 2)
                if not self.useDist:
                    sim.checkpoint(hdfprefix)
                    sim.solver.reset()
                    sim.restore(hdfprefix)
                sim.run(self.endTime)

        if MPI._shouldWrite:
            with HDF5Handler(hdfprefix) as hdf:
                res, = hdf['test'].results

                minVals = np.min(res.data, axis=0)
                maxVals = np.max(res.data, axis=0)

                # Plotting code in case manual inspection is needed
                if plot:
                    from matplotlib import pyplot as plt
                    ntot = res.data[...].shape[2]
                    nrows = int(ntot ** 0.5)
                    ncols = int(np.ceil(ntot / nrows))
                    for i in range(res.data[...].shape[2]):
                        plt.subplot(nrows, ncols, i + 1)
                        plt.plot(res.time[0], res.data[:,:,i].T)
                        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                        plt.title(res.labels[i])
                    plt.show()

                # Check that all runs have exactly the same traces
                self.assertTrue((minVals == maxVals).all())


class TetexactReproducibilityTestCase(ReproducibilityTestCase, base_model.TetTestModelFramework):
    """Test reproducibility of Tetexact"""

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation('Tetexact', nmdl, ngeom, nrng, True)
        return nsim

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testReproducible(self):
        self._testReproducible(plot=False)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TetexactReproducibilityTestCase))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())


