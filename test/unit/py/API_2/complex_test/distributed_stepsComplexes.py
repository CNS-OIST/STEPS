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

"""
Tests for steps complex reaction simulation with distTetOpSplit. Check that steps complex reactions
give the same results as complexes declared with statesAsSpecies=True.
"""

import numpy as np
import os
import operator
import unittest

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

from . import test_stepsComplexes

TEST_DIR = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '..', '..', '..', '..')
MESH_DIR = os.path.join(TEST_DIR, 'mesh')


class StepsDistComplexReaction(test_stepsComplexes.StepsComplexReaction):
    """Test Complex reactions with distTetOpSplit."""

    def setUp(self, *args, **kwargs):
        self.membs = []
        super().setUp(*args, **kwargs, complexCls=Channel)

    def getSimulation(self, mdl, geom, rng):
        # TODO: Activate once complex diffusion is implemented
        # with mdl, mdl.vsys:
        #     if mdl.Comp1._areStatesAsSpecies():
        #         Diffusion(mdl.Comp1, 1e-15)  # 1e-14)  #
        #         Diffusion(mdl.Comp2, 1e-15)  # 1e-14)  #
        return Simulation('DistTetOpSplit', mdl, geom, rng)

    def getGeometries(self):
        meshPath = os.path.join(MESH_DIR, 'cube_379tets.msh')
        geoms = [DistMesh(meshPath, 0.75e-6), DistMesh(meshPath, 0.75e-6)]
        for geom, vsys, ssys in zip(geoms, self.volsyss, self.surfsyss):
            # Construct compartment so that patch triangles are all owned by the same process
            with geom.asLocal():
                if MPI.rank == 0:
                    localTets = TetList(geom.tets)
                    for tri in localTets.surface - geom.surface:
                        localTets -= tri.tetNeighbs
                else:
                    localTets = TetList(local=True)
                compTets = localTets.combineWithOperator(operator.or_)

            with geom:
                cyt = Compartment.Create(compTets, vsys, conductivity=1)
                out = Compartment.Create(geom.tets - cyt.tets, vsys, conductivity=1)
                self.comps.append(cyt)
                self.comps.append(out)

                patch = Patch.Create(cyt.surface & out.surface, cyt, out, ssys)
                self.patches.append(patch)

                memb = Membrane.Create([patch])
                self.membs.append(memb)
        return geoms

    def initializeSimulation(self, sim):
        super().initializeSimulation(sim)
        sim.memb.Potential = -0.07
        sim.EfieldDT = 1e-1
        sim.TRIS(sim.geom.patch.tris).Capac = 0.01
        
    def _testVDepComplexReacs(self, same):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, self.surfsyss, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    v1 = VDepRate(lambda V: 1 / (1 + np.exp(-(-0.07 - V) / 0.01)))
                    v2 = VDepRate(lambda V: 1 / (1 + np.exp((-0.04 - V) / 0.01)))
                    with Comp1[..., sus2A, sus2A]:
                        sus1A.s <r[1]> sus1B.s <r[2]> sus1C.s
                        if same or mdl is self.mdls[0]:
                            r[1].K = v1, v2
                            r[2].K = 2 * v1, 2 * v2
                        else:
                            r[1].K = 5 * v1, 2 * v2
                            r[2].K = 10 * v1, 3 * v2
                    curr1 = OhmicCurr.Create(Comp1[sus1C, ...], 1e-16, 0.02)
                    curr2 = OhmicCurr.Create(Comp1[sus1A, sus1A, ...], 1e-15, -0.07)

            sim = self.getSimulation(mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.patch.LIST(*Comp1[..., sus2A, sus2A]).Count
            data <<= rs.patch.Comp1[..., sus2A, sus2A].sus1A.Count
            data <<= rs.SUM(rs.TRIS(geom.patch.tris).V) / len(geom.patch.tris)
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                self.initializeSimulation(sim)
                sim.patch.Comp1[sus1A, sus1A, sus2A, sus2A].Count = 250
                sim.run(self.ENDT // 2)
                sim.memb.Potential = -0.03
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testVDepComplexReacs(self):
        self._testVDepComplexReacs(same=True)
        self.setUp()
        self._testVDepComplexReacs(same=False)

