####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

""" Unit tests for the parameter setting and export system."""

import filecmp
import os
import tempfile
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.utils import *

from . import base_model

VALID_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parameter_export')


class TetParamExport(base_model.TetTestModelFramework):
    """Test parameter declaration."""
    def setUp(self):
        super().setUp()

        self._initC1S1 = self.initC1S1
        self._r01f = self.r01f
        self._r02f = self.r02f
        self._Chan1_G = self.Chan1_G
        self._Chan1_rev = self.Chan1_rev

        # Declare Parameters
        self.initC1S1 = Parameter(self._initC1S1, name='initC1S1', DOI='test DOI', otherComment='comment with spaces')
        self.initC1CCsus11 = Parameter(self.initC1CCsus11, name='initC1CCsus11', DOI2='other DOI')
        self.r01f = Parameter(self._r01f * 1e-6, 'uM^-1 s^-1', name='r01f')
        self.r02f = Parameter(self._r02f / 2 * 1e-6, 'uM^-1 s^-1', name='r02f') * 2 + self.r01f - self.r01f
        self.initP1Ex = Parameter(self.initP1Ex, name='initP1Ex')
        self.Chan1_G = Parameter(self._Chan1_G * 1e6, 'uS', name='Chan1_G')
        self.Chan1_rev = Parameter(self._Chan1_rev * 1e3, 'mV', name='Chan1_rev')

        self.mdl = self.get_API2_Mdl()
        self.mesh = self.get_API2_Geom(self.mdl)

        self.sim = self._get_API2_Sim(self.mdl, self.mesh)

        self.createdFiles = []

    def tearDown(self):
        super().tearDown()
        for path in self.createdFiles:
            if os.path.isfile(path):
                os.remove(path)

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        nsim = Simulation.Create('Tetexact', nmdl, ngeom, nrng, True)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _testExport(self, testName, **kwargs):
        _, filePath = tempfile.mkstemp(prefix=testName)
        self.createdFiles.append(filePath)

        # Do not check the Function column as it shows full paths to python files
        exportedPaths = ExportParameters(self.sim, filePath, hideColumns=['Function'], **kwargs)
        self.createdFiles += exportedPaths

        validPaths = [
            fname[len(testName):] for fname in os.listdir(os.path.join(VALID_DIR)) if fname.startswith(testName)
        ]

        self.assertEqual(
            set(validPaths),
            set([os.path.basename(path)[len(os.path.basename(filePath)):] for path in exportedPaths]),
        )

        for path in exportedPaths:
            specificName = path[len(filePath):]
            validPath = os.path.join(VALID_DIR, testName + specificName)

            self.assertTrue(filecmp.cmp(path, validPath, shallow=False))

    def testDefaultExport(self):
        self._testExport(f'{self.__class__.__name__}_defaultExport')

    def testDefaultLatexExport(self):
        self._testExport(f'{self.__class__.__name__}_defaultLatexExport', textFormat='latex')

    def testInitExport(self):
        self.sim.newRun()
        self.init_API2_sim(self.sim)
        self._testExport(f'{self.__class__.__name__}_initExport')

    def testParamSetting(self):
        self.sim.newRun()
        self.init_API2_sim(self.sim)

        self.assertEqual(
            self.sim.comp1.S1.Count,
            self._initC1S1,
        )
        self.assertEqual(
            self.mdl.vs1R1['fwd']._getStepsObjects()[0].getKcst(),
            self._r01f,
        )
        self.assertEqual(
            self.mdl.vs2R2['fwd']._getStepsObjects()[0].getKcst(),
            self._r02f
        )
        self.assertEqual(
            self.mdl.ssys.Chan1_Ohm_I._getStepsObjects()[0].getG(),
            self._Chan1_G
        )
        self.assertEqual(
            self.mdl.ssys.Chan1_Ohm_I._getStepsObjects()[0].getERev(),
            self._Chan1_rev
        )


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TetParamExport, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

