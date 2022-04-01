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

""" Unit tests for Tetexact small current recording."""

import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.rng import *
from steps.saving import *

import os

FILEDIR = os.path.dirname(os.path.abspath(__file__))

OUT_CONC = 2e-3
CHAN_P = 2.5e-19

POT = -60e-3
RA = 1.2
CAP = 0.64e-2

SIM_END = 1e-2
EF_DT = 1e-7
INT_THRESH = 0.001

QE = 1.602176634e-19

class SmallCurrentSavingTest(unittest.TestCase):
    """Test whether the bounding box is computed correctly with split meshes."""
    def getSimulation(self):
        mdl = Model()
        with mdl:
            ssys = SurfaceSystem.Create()
            ion = Species.Create(valence=1)
            
            ss1 = SubUnitState.Create()
            chan = Channel.Create([ss1])
            
            with ssys:
                curr = GHKCurr.Create(chan[ss1], ion, CHAN_P, virtual_oconc=OUT_CONC, computeflux=True)
            
        mesh = TetMesh.LoadGmsh(os.path.join(FILEDIR, 'meshes', '2_tets.msh'), 1e-6)
        with mesh:
            comp = Compartment.Create(mesh.tets[:1])
            tris = comp.surface & (mesh.tets - comp.tets).surface
            patch = Patch.Create(tris, comp, None, ssys)
            memb = Membrane.Create([patch])

        rng = RNG('mt19937', 512, 132)
        return Simulation('Tetexact', mdl, mesh, rng, True), tris

    def runSimulation(self, sim):
        sim.newRun()
        sim.EfieldDT = EF_DT
        sim.memb.Potential = POT
        sim.memb.VolRes = RA
        sim.memb.Capac = CAP

        sim.patch.chan[sim.model.ss1].Count = 1

        sim.run(SIM_END)

    def getNbCorrectCurrent(self, res):
        correct = 0
        # Check whether currents correspond to an integer number of ions
        for val in res.data[0, :, 0]:
            qv = val / QE * EF_DT
            if abs(qv - round(qv)) < INT_THRESH:
                correct += 1
        return correct

    def testSimdtEqualEfdt(self):
        sim, tris = self.getSimulation()
        rs = ResultSelector(sim)

        res = rs.TRIS(tris).curr.I << rs.comp.ion.Count
        sim.toSave(res, dt=EF_DT)

        self.runSimulation(sim)

        self.assertEqual(self.getNbCorrectCurrent(res), len(res.time[0]))

        # Check that the current corresponds to each ion change
        for i in range(len(res.time[0])):
            qv = res.data[0, i, 0] / QE * EF_DT
            dn = res.data[0, i, 1] - (res.data[0, i - 1, 1] if i > 0 else 0)
            self.assertLess(abs(qv + dn), INT_THRESH)

    def testSimdtEqual2Efdt(self):
        sim, tris = self.getSimulation()
        rs = ResultSelector(sim)

        res = rs.TRIS(tris).curr.I << rs.comp.ion.Count
        sim.toSave(res, dt=2 * EF_DT)

        self.runSimulation(sim)

        self.assertEqual(self.getNbCorrectCurrent(res), len(res.time[0]))

        # Check that the current is lower than the maximum current for each ion change
        for i in range(len(res.time[0])):
            qv = res.data[0, i, 0] / QE * EF_DT
            dn = res.data[0, i, 1] - (res.data[0, i - 1, 1] if i > 0 else 0)
            if dn < 0:
                self.assertLess(dn - INT_THRESH, qv)
            else:
                self.assertGreater(dn + INT_THRESH, -qv)

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(SmallCurrentSavingTest, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
