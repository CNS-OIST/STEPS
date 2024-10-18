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

"""
Tests for steps complex reaction simulation. Check that steps complex reactions
give the same results as complexes declared with statesAsSpecies=True.
"""

import unittest

import scipy.stats
import numpy as np

from steps import interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

import time

class StepsComplexReaction(unittest.TestCase):
    """Test Complex creation."""
    def setUp(self, order=NoOrdering, nbRuns=50, pThresh=0.001):
        self.nbRuns = nbRuns
        self.pThresh = pThresh
        self.ENDT = 5

        self.mdls = [Model(), Model()]
        self.geoms = [Geometry(), Geometry()]
        self.r = ReactionManager()

        self.sass = [True, False]
        self.species = []
        self.volsyss = []
        self.surfsyss = []
        self.complexes = []
        self.suss = []
        self.subs = []
        self.comps = []
        self.patches = []
        for mdl, sas in zip(self.mdls, self.sass):
            with mdl:
                SA, SB = Species.Create()
                self.species.append([SA, SB])

                sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = SubUnitState.Create()
                self.suss.append([sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C])

                su1, su2, su3, su4 = SubUnit.Create(
                    [sus1A, sus1B, sus1C], 
                    [sus2A, sus2B, sus2C],
                    [sus3A, sus3B, sus3C],
                    [sus4A, sus4B, sus4C],
                )
                self.subs.append([su1, su2, su3, su4])

                Comp1 = Complex.Create([su1]*2 + [su2]*2, statesAsSpecies=sas, order=order)
                Comp2 = Complex.Create([su3]*2 + [su4]*2, statesAsSpecies=sas, order=order)
                self.complexes.append([Comp1, Comp2])

                vsys = VolumeSystem.Create()
                self.volsyss.append(vsys)

                ssys = SurfaceSystem.Create()
                self.surfsyss.append(ssys)

        for geom, vsys, ssys in zip(self.geoms, self.volsyss, self.surfsyss):
            with geom:
                cyt = Compartment.Create(vsys, 1.6572e-19)
                out = Compartment.Create(vsys, 1.6572e-19)

                self.comps.append(cyt)
                self.comps.append(out)

                patch = Patch.Create(cyt, out, ssys, 1e-12)
                self.patches.append(patch)

        # Hardcoded seed to avoid random failures
        self.rng = RNG('mt19937', 512, 123456)

    @staticmethod
    def _addSaverToSim(sim, locations, complexes, addSingleSUS=True):
        rs = ResultSelector(sim)
        data = rs.LIST(*locations).ALL(Species).Count
        for comp in complexes:
            data <<= rs.LIST(*locations).LIST(*comp).Count
            if addSingleSUS:
                data <<= rs.LIST(*locations).LIST(comp).LIST(list(comp._subUnits[0]._states)[0]).Count
        sim.toSave(data, dt=0.1)
        return data

    def _testData(self, allData, same=True, labels=None, plot=False):
        # Perform Kolmogorov-Smirnov two sampled tests on each time point
        dat1, dat2 = allData
        self.assertEqual(dat1.shape, dat2.shape)
        nTot = dat1.shape[2]

        if plot:
            from matplotlib import pyplot as plt
            nRows = int(nTot ** 0.5)
            nCols = int(np.ceil(nTot / nRows))
            fig = None

        totNbFail = 0
        for j in range(nTot):
            nbTests = 0
            allPVals = []
            for t in range(dat1.shape[1]):
                d, p = scipy.stats.mstats.ks_twosamp(dat1[:,t,j], dat2[:,t,j])
                nbTests += 1
                allPVals.append(p)

            # Raise an error if the number of "failed" tests is higher than expected
            # Apply Benjamini–Hochberg procedure to get the number of "true discoveries"
            nbFail = nbTests
            for k, pval in enumerate(sorted(allPVals)):
                if pval > (k + 1) / nbTests * self.pThresh:
                    nbFail = k
                    break

            totNbFail += min(nbFail, 1)

            # Plot average time traces of values that are significantly different from expected
            if plot and ((same and nbFail > 0) or (not same and nbFail == 0)):
                if fig is None:
                    fig = plt.figure()
                    fig.suptitle(f'Data should be ' + ('identical' if same else 'different'))
                plt.subplot(nRows, nCols, j+1)
                plt.plot(range(dat1.shape[1]), np.mean(dat1[:,:,j],axis=0), 'b', range(dat1.shape[1]), np.mean(dat2[:,:,j],axis=0), 'r')
                plt.plot(range(dat1.shape[1]), np.log(allPVals), 'g')
                if labels is not None:
                    plt.title(labels[j])

        if plot and fig is not None:
            if (same and totNbFail > 0) or (not same and totNbFail == 0):
                plt.show()
            else:
                plt.close()

        if same:
            self.assertEqual(totNbFail, 0)
        else:
            self.assertGreater(totNbFail, 0)

    def testCreation(self):
        mdl = Model()
        with mdl:
            sus1, sus2 = SubUnitState.Create()
            su = SubUnit.Create([sus1, sus2])
            with self.assertRaises(Exception):
                CC = Complex.Create([su], statesAsSpecies=False, order=StrongOrdering)
    
    def _runComplexSurfaceReacSubTests(self, func):
        with self.subTest('Volume-Volume reactions'):
            func(self.volsyss, lambda x, _: x, ['cyt', 'cyt'], True)
            self.setUp()
            func(self.volsyss, lambda x, _: x, ['cyt', 'cyt'], False)
        with self.subTest('Surface-Surface reactions'):
            self.setUp()
            func(self.surfsyss, lambda x, _: Surf(x), ['patch', 'patch'], True, kfact=2e3)
            self.setUp()
            func(self.surfsyss, lambda x, _: Surf(x), ['patch', 'patch'], False, kfact=2e3)
        with self.subTest('Volume-Surface reactions'):
            self.setUp()
            func(self.surfsyss, lambda x, i: [Surf, In][i](x), ['patch', 'cyt'], True)
            self.setUp()
            func(self.surfsyss, lambda x, i: [Surf, In][i](x), ['patch', 'cyt'], False)
            self.setUp()
            func(self.surfsyss, lambda x, i: [In, Surf][i](x), ['cyt', 'patch'], True)
            self.setUp()
            func(self.surfsyss, lambda x, i: [In, Surf][i](x), ['cyt', 'patch'], False)

    def _testSubUnitChange(self, systems, locFunc, locations, same=True, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    with Comp1[..., sus2A, sus2A]:
                        locFunc(sus1B, 0) <r[1]> locFunc(sus1A, 0) <r[2]> locFunc(sus1C, 0)
                        locFunc(sus1A['a'], 0) + locFunc(sus1B['b'], 0) <r[3]> locFunc(sus1C['a'], 0) + locFunc(sus1C['b'], 0)
                        if same or mdl is self.mdls[0]:
                            r[1].K = 1, 1
                            r[2].K = 1, 1
                            r[3].K = 2, 2
                        else:
                            r[1].K = 5, 1
                            r[2].K = 2, 1
                            r[3].K = 3, 2

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[0]).LIST(*Comp1[..., sus2A, sus2A]).Count
            data <<= rs.LIST(locations[0]).Comp1[..., sus2A, sus2A].sus1A.Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2A, sus2A].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testSubUnitChange(self):
        self._runComplexSurfaceReacSubTests(self._testSubUnitChange)

    def _testSubUnitSelectors(self, systems, locFunc, locations, same=True, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    with Comp1[..., sus2A, sus2A]:
                        locFunc(sus1B | sus1A, 0) >r[1]> locFunc(sus1C, 0)
                        locFunc(sus1C | sus1A, 0) >r[2]> locFunc(sus1B, 0)
                        locFunc(sus1C | sus1B, 0) >r[3]> locFunc(sus1A, 0)
                        locFunc(sus1A | sus1B, 0) + locFunc(SA, 1) >r[4]> locFunc(sus1C, 0)
                        if same or mdl is self.mdls[0]:
                            r[1].K = 5
                            r[2].K = 3
                            r[3].K = 1
                            r[4].K = 1000
                        else:
                            r[1].K = 10
                            r[2].K = 1
                            r[3].K = 2
                            r[4].K = 2000

                    # Check that subunitstates in filters are not counted as reactants
                    with Comp1[sus1A, ..., sus2A, sus2A]:
                        locFunc(sus1A, 0) <r[1]> locFunc(sus1C, 0)
                        if same or mdl is self.mdls[0]:
                            r[1].K = 2, 3
                        else:
                            r[1].K = 1, 5

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[1]).SA.Count
            data <<= rs.LIST(locations[0]).LIST(*Comp1[..., sus2A, sus2A]).Count
            data <<= rs.LIST(locations[0]).Comp1[..., sus2A, sus2A].sus1A.Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[1]).SA.Count = 400
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2A, sus2A].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testSubUnitSelectors(self):
        self._runComplexSurfaceReacSubTests(self._testSubUnitSelectors)

    def _testSpeciesInteraction(self, systems, locFunc, locations, same=True, kfact=1, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    with Comp1[..., sus2A, sus2A]:
                        locFunc(sus1A, 0) + locFunc(SA, 1) <r[1]> locFunc(sus1B, 0) + locFunc(SA, 1) <r[2]> locFunc(sus1C, 0) + locFunc(SB, 1)
                        locFunc(sus1C, 0) + locFunc(SB, 0) <r[3]> locFunc(sus1B, 1) + locFunc(SB, 0)
                        if same or mdl is self.mdls[0]:
                            r[1].K = kfact*10000, kfact*1000
                            r[2].K = kfact*20000, kfact*2000
                            r[3].K = kfact*5000, kfact*1000
                        else:
                            r[1].K = kfact*2000, kfact*7000
                            r[2].K = kfact*40000, kfact*12000
                            r[3].K = kfact*20000, kfact*7000
            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[1]).LIST(SA, SB).Count
            data <<= rs.LIST(locations[0]).LIST(*Comp1[..., sus2A, sus2A]).Count
            if locations[0] != locations[1]:
                data <<= rs.LIST(locations[1]).LIST(Comp1[sus1B, sus1B, sus2A, sus2A]).Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[1]).SA.Count = 400
                sim.LIST(locations[1]).SB.Count = 400
                sim.LIST(locations[0]).SB.Count = 400
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2A, sus2A].Count = 500
                sim.LIST(locations[1]).Comp1[sus1B, sus1B, sus2A, sus2A].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testSpeciesInteraction(self):
        self._runComplexSurfaceReacSubTests(self._testSpeciesInteraction)

    def _testTwoComplexes(self, systems, locFunc, locations, same=True, kfact=1, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    with Comp1[..., sus2A, sus2A], Comp2[..., sus4A, sus4A]:
                        locFunc(sus1A, 0) + locFunc(sus3A, 1) <r[1]> locFunc(sus1A, 0) + locFunc(sus3B, 1)
                        locFunc(sus1A, 0) + locFunc(sus3B, 1) <r[2]> locFunc(sus1B, 0) + locFunc(sus3B, 1)
                        locFunc(sus1B, 0) + locFunc(sus3B, 1) <r[3]> locFunc(sus1B, 0) + locFunc(sus3C, 1)
                        locFunc(sus1B, 0) + locFunc(sus3C, 1) <r[4]> locFunc(sus1C, 0) + locFunc(sus3C, 1)
                        if same or mdl is self.mdls[0]:
                            r[1].K = kfact * 10000, kfact * 10000
                            r[2].K = kfact * 10000, kfact * 10000
                            r[3].K = kfact * 10000, kfact * 10000
                            r[4].K = kfact * 10000, kfact * 10000
                        else:
                            r[1].K = kfact * 20000, kfact * 10000
                            r[2].K = kfact * 5000, kfact * 10000
                            r[3].K = kfact * 15000, kfact * 10000
                            r[4].K = kfact * 1000, kfact * 10000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[0]).LIST(*Comp1[..., sus2A, sus2A]).Count
            data <<= rs.LIST(locations[1]).LIST(*Comp2[..., sus4A, sus4A]).Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2A, sus2A].Count = 500
                sim.LIST(locations[1]).Comp2[sus3A, sus3A, sus4A, sus4A].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testTwoComplexes(self):
        self._runComplexSurfaceReacSubTests(self._testTwoComplexes)

    def _testTwoSameComplexes(self, systems, locFunc, locations, same=True, kfact=1, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    with Comp1[..., sus2C, sus2C] as C1, Comp1[..., sus2C, sus2C] as C2:
                        locFunc(sus1A[C1], 0) + locFunc(sus1A[C2], 1) <r[1]> locFunc(sus1A[C1], 0) + locFunc(sus1B[C2], 1)
                        locFunc(sus1A[C1], 0) + locFunc(sus1B[C2], 1) <r[2]> locFunc(sus1B[C1], 0) + locFunc(sus1B[C2], 1)
                        locFunc(sus1B[C1], 0) + locFunc(sus1B[C2], 1) <r[3]> locFunc(sus1B[C1], 0) + locFunc(sus1C[C2], 1)
                        locFunc(sus1B[C1], 0) + locFunc(sus1C[C2], 1) <r[4]> locFunc(sus1C[C1], 0) + locFunc(sus1C[C2], 1)
                        if same or mdl is self.mdls[0]:
                            r[1].K = kfact * 10000, kfact * 10000
                            r[2].K = kfact * 10000, kfact * 10000
                            r[3].K = kfact * 10000, kfact * 10000
                            r[4].K = kfact * 10000, kfact * 10000
                        else:
                            r[1].K = kfact * 20000, kfact * 10000
                            r[2].K = kfact * 5000, kfact * 10000
                            r[3].K = kfact * 15000, kfact * 10000
                            r[4].K = kfact * 1000, kfact * 10000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[0]).LIST(*Comp1[..., sus2C, sus2C]).Count
            if locations[0] != locations[1]:
                data <<= rs.LIST(locations[1]).LIST(*Comp1[..., sus2C, sus2C]).Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2C, sus2C].Count = 500
                sim.LIST(locations[1]).Comp1[sus1A, sus1A, sus2C, sus2C].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testTwoSameComplexes(self):
        self._runComplexSurfaceReacSubTests(self._testTwoSameComplexes)

    def _testFullcomplexes(self, systems, locFunc, locations, same=True, kfact=1, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    C1, C2 = Comp1.get(), Comp2.get()

                    locFunc(C1[sus1A, :, sus2C, sus2C], 0) + locFunc(C2[sus3A, :, sus4C, sus4C], 1) <r[1]> locFunc(C1[sus1A, :, sus2C, sus2C], 0) + locFunc(C2[sus3B, :, sus4C, sus4C], 1)
                    if same or mdl is self.mdls[0]:
                        r[1].K = kfact * 10000, kfact * 10000
                    else:
                        r[1].K = kfact * 20000, kfact * 10000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[0]).LIST(*Comp1[sus1A, :, sus2C, sus2C]).Count
            data <<= rs.LIST(locations[1]).LIST(*Comp2[~sus3C, :, sus4C, sus4C]).Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2C, sus2C].Count = 500
                sim.LIST(locations[1]).Comp2[sus3A, sus3A, sus4C, sus4C].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testFullcomplexes(self):
        self._runComplexSurfaceReacSubTests(self._testFullcomplexes)

    def _testFullcomplexesSubSelectors(self, systems, locFunc, locations, same=True, kfact=1, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    C1 = Comp1.get()

                    locFunc(C1[sus1A|sus1B, :, sus2C, sus2C], 0) >r[1]> locFunc(C1[sus1C, :, sus2C, sus2C], 0)
                    locFunc(C1[sus1A|sus1C, :, sus2C, sus2C], 0) >r[2]> locFunc(C1[sus1B, :, sus2C, sus2C], 0)
                    locFunc(C1[~sus1A, :, sus2C, sus2C], 0) >r[3]> locFunc(C1[sus1A, :, sus2C, sus2C], 1)
                    locFunc(C1[~sus1B, :, sus2C, sus2C], 0) >r[4]> locFunc(C1[sus1B, :, sus2C, sus2C], 1)
                    if same or mdl is self.mdls[0]:
                        r[1].K = 4
                        r[2].K = 3
                        r[3].K = 2
                        r[4].K = 1
                    else:
                        r[1].K = 1
                        r[2].K = 2
                        r[3].K = 3
                        r[4].K = 4

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[0]).LIST(*Comp1[..., sus2C, sus2C]).Count
            if locations[0] != locations[1]:
                data <<= rs.LIST(locations[1]).LIST(*Comp1[~sus1C, :, sus2C, sus2C]).Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[0]).Comp1[sus1A, sus1B, sus2C, sus2C].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testFullcomplexesSubSelectors(self):
        self._runComplexSurfaceReacSubTests(self._testFullcomplexesSubSelectors)

    def _testComplexCreationDeletion(self, systems, locFunc, locations, same=True, kfact=1, **kwargs):
        r = self.r
        allData = []
        for mdl, geom, sys, Comp, subunitstates, species in zip(self.mdls, self.geoms, systems, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with sys:
                    with Comp1[...] as C1:
                        locFunc(sus1A, 0) + locFunc(SA, 1) <r[1]> locFunc(sus1B, 0) + locFunc(SA, 1)
                        locFunc(sus1B, 0) + locFunc(SB, 1) <r[2]> locFunc(sus1C, 0) + locFunc(SB, 1)
                        if same or mdl is self.mdls[0]:
                            r[1].K = kfact * 80000, kfact * 10000
                            r[2].K = kfact * 50000, kfact * 10000
                        else:
                            r[1].K = kfact * 10000, kfact * 50000
                            r[2].K = kfact * 30000, kfact * 10000

                    C1 = Comp1.get()

                    locFunc(C1[sus1B, sus1B, ...], 0) >r[1]> None
                    if same or mdl is self.mdls[0]:
                        r[1].K = 10
                    else:
                        r[1].K = 5

                    locFunc(C1[sus1C, sus1C, ...], 0) + locFunc(SB, 0) >r[1]> locFunc(C1[sus1C, sus1C, ...], 0) + locFunc(SB, 0) + locFunc(Comp1[sus1A, sus1A, sus2A, sus2A], 0)
                    if same or mdl is self.mdls[0]:
                        r[1].K = kfact * 1000000
                    else:
                        r[1].K = kfact * 1500000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            rs = ResultSelector(sim)
            data = rs.LIST(locations[0]).LIST(*Comp1[..., sus2C, sus2C]).Count
            sim.toSave(data, dt=0.1)

            for i in range(self.nbRuns):
                sim.newRun()
                sim.LIST(locations[1]).SA.Count = 500
                sim.LIST(locations[1]).SB.Count = 500
                sim.LIST(locations[0]).Comp1[sus1A, sus1A, sus2C, sus2C].Count = 500
                sim.run(self.ENDT)

            allData.append(np.array(data.data[:]))

        self._testData(allData, same=same, labels=data.labels)

    def testComplexCreationDeletion(self):
        self._runComplexSurfaceReacSubTests(self._testComplexCreationDeletion)

    def testComplexDatafromSim(self):
        r = self.r
        allData = []
        for mdl, geom, vsys, Comp, subunitstates, species in zip(self.mdls, self.geoms, self.volsyss, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with vsys:
                    SA <r[1]> SB
                    r[1].K = 1, 1

                    with Comp1[...] as C1:
                        sus1A + SA <r[1]> sus1B + SA <r[2]> sus1C + SA
                        r[1].K = 20000, 10000
                        r[2].K = 10000, 20000

                    with Comp2[...] as C2:
                        sus3A + SA <r[1]> sus3B + SA <r[2]> sus3C + SA
                        r[1].K = 20000, 10000
                        r[2].K = 10000, 20000


            sim = Simulation('Wmdirect', mdl, geom, self.rng)

            sim.cyt.Comp1.sus1A.Count
            with self.assertRaises(Exception):
                sim.cyt.Comp1.sus3A.Count


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(StepsComplexReaction))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
