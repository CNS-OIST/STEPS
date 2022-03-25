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
    #def setUp(self, order=NoOrdering, nbRuns=100, pThresh=0.01):
    def setUp(self, order=NoOrdering, nbRuns=100, pThresh=0.001):
        self.nbRuns = nbRuns
        self.pThresh = pThresh

        self.mdls = [Model(), Model()]
        self.geoms = [Geometry(), Geometry()]
        self.r = ReactionManager()

        self.sass = [True, False]
        self.species = []
        self.volsyss = []
        self.complexes = []
        self.suss = []
        self.subs = []
        self.comps = []
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

        for geom, vsys in zip(self.geoms, self.volsyss):
            with geom:
                # Create the cytosol and Endoplasmic Reticulum compartments
                cyt = Compartment.Create(vsys, 1.6572e-19)

                self.comps.append(cyt)

        self.rng = RNG('mt19937', 512, time.time())

    @staticmethod
    def _addSaverToSim(sim, comps, addSingleSUS=True):
        rs = ResultSelector(sim)
        for comp in comps:
            states = list(comp)
            comp = comp._getReferenceObject()
            data = getattr(rs.cyt, comp.name)[tuple(states[0]._state)].Count
            for s in states[1:]:
                data <<= getattr(rs.cyt, comp.name)[tuple(s._state)].Count
            if addSingleSUS:
                data <<= getattr(getattr(rs.cyt, comp.name), states[0]._state[0].name).Count
        data <<= rs.cyt.SA.Count
        data <<= rs.cyt.SB.Count
        sim.toSave(data, dt=0.1)
        return data

    def _testData(self, allData, same=True):
        # Perform Kolmogorov-Smirnov two sampled tests on each time point
        dat1, dat2 = allData
        self.assertEqual(dat1.shape, dat2.shape)
        totNbFail = 0
        for j in range(dat1.shape[2]):
            nbTests = 0
            allPVals = []
            for t in range(dat1.shape[1]):
                d, p = scipy.stats.mstats.ks_twosamp(dat1[:,t,j], dat2[:,t,j])
                nbTests += 1
                allPVals.append(p)

            # Raise an error if the number of "failed" tests is higher than expected
            # Apply Benjamini–Hochberg procedure to get the number of "true discoveries"
            allPVals.sort()
            nBFail = nbTests
            for k, pval in enumerate(allPVals):
                if pval > (k + 1) / nbTests * self.pThresh:
                    nbFail = k
                    break

            totNbFail += min(nbFail, 1)

            # Import pylab and uncomment to plot average time traces of values that are significantly different

            # if same and nbFail > 0:
                # print('fail', j)
                # pylab.plot(range(dat1.shape[1]), np.mean(dat1[:,:,j],axis=0), 'b', range(dat1.shape[1]), np.mean(dat2[:,:,j],axis=0), 'r')
                # pylab.show()

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

    def testSubUnitChange(self):
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

                    with Comp1[...]:
                        sus1B <r[1]> sus1A <r[2]> sus1C
                        sus1A['a'] + sus1B['b'] <r[3]> sus1C['a'] + sus1C['b']
                        r[1].K = 1, 1
                        r[2].K = 1, 1
                        r[3].K = 2, 2

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testSubUnitSelectors(self):
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

                    with Comp1[...]:
                        (sus1B | sus1A) >r[1]> sus1C
                        (sus1C | sus1A) >r[2]> sus1B
                        (sus1C | sus1B) >r[3]> sus1A
                        r[1].K = 5
                        r[2].K = 3
                        r[3].K = 1

                    # Check that subunitstates in filters are not counted as reactants
                    with Comp1[sus1A, ...]:
                        sus1A <r[1]> sus1C
                        r[1].K = 2, 3

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testSubUnitChangeDiffParams(self):
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

                    with Comp1[...]:
                        sus1B <r[1]> sus1A <r[2]> sus1C
                        if mdl is self.mdls[0]:
                            r[1].K = 1, 1
                            r[2].K = 1, 1
                        else:
                            r[1].K = 100.5, 1.1
                            r[2].K = 0.7, 1.2

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData, same=False)

    def testSpeciesInteraction(self):
        r = self.r
        allData = []
        for mdl, geom, vsys, Comp, subunitstates, species in zip(self.mdls, self.geoms, self.volsyss, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            with mdl:
                with vsys:
                    with Comp1[...]:
                        sus1A + SA <r[1]> sus1B + SA <r[2]> sus1C + SB
                        if mdl is self.mdls[0]:
                            r[1].K = 10000, 15000
                        else:
                            r[1].K = 10000, 15000
                        r[2].K = 20000, 25000
            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testTwoComplexes(self):
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

                    with Comp1[...], Comp2[...]:
                        sus1A + sus3A <r[1]> sus1A + sus3B <r[2]> sus1B + sus3B <r[3]> sus1B + sus3C <r[4]> sus1C + sus3C
                        r[1].K = 10000, 10000
                        r[2].K = 10000, 10000
                        r[3].K = 10000, 10000
                        r[4].K = 10000, 10000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1, Comp2])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.cyt.Comp2[(sus3A,)*2 + (sus4C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testTwoSameComplexes(self):
        r = self.r
        allData = []
        for mdl, geom, vsys, Comp, subunitstates, species in zip(self.mdls, self.geoms, self.volsyss, self.complexes, self.suss, self.species):
            Comp1, Comp2 = Comp
            sus1A, sus1B, sus1C, sus2A, sus2B, sus2C, sus3A, sus3B, sus3C, sus4A, sus4B, sus4C = subunitstates
            SA, SB = species
            #print('Run with sas = ', Comp1._statesAsSpecies)
            with mdl:
                with vsys:
                    SA <r[1]> SB
                    r[1].K = 1, 1

                    with Comp1[..., sus2C, sus2C] as C1, Comp1[..., sus2C, sus2C] as C2:
                        sus1A[C1] + sus1A[C2] <r[1]> sus1A[C1] + sus1B[C2] <r[2]> sus1B[C1] + sus1B[C2] <r[3]> sus1B[C1] + sus1C[C2] <r[4]> sus1C[C1] + sus1C[C2]
                        r[1].K = 10000, 10000
                        r[2].K = 10000, 10000
                        r[3].K = 10000, 10000
                        r[4].K = 10000, 10000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1[..., sus2C, sus2C]], addSingleSUS=False)

            for i in range(self.nbRuns):
                sim.newRun()
                #print('newRun', i)
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                #print('start run')
                sim.run(10)
                #print('end run')

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testFullcomplexes(self):
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

                    C1, C2 = Comp1.get(), Comp2.get()

                    C1[...] <r[1]> C1[...]
                    C2[...] <r[2]> C2[...]
                    r[1].K = 1, 1
                    r[2].K = 1, 1

                    C1[sus1A, ..., sus2C] + C2[sus3A, ..., sus4C] <r[1]> C1[sus1A, ..., sus2C] + C2[sus3B, ..., sus4C]
                    r[1].K = 10000, 10000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1, Comp2])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.cyt.Comp2[(sus3A,)*2 + (sus4C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testFullcomplexesSubSelectors(self):
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

                    C1 = Comp1.get()

                    C1[...] <r[1]> C1[...]
                    r[1].K = 1, 1

                    C1[sus1A|sus1B, ..., sus2C] >r[1]> C1[sus1C, ..., sus2C]
                    C1[sus1A|sus1C, ..., sus2C] >r[2]> C1[sus1B, ..., sus2C]
                    r[1].K = 1
                    r[2].K = 1

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

    def testComplexCreationDeletion(self):
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

                    C1 = Comp1.get()

                    C1[sus1B, sus1B, ...] >r[1]> None
                    r[1].K = 10

                    C1[sus1C, sus1C, ...] + SB >r[1]> C1[sus1C, sus1C, ...] + SB + Comp1[sus1A, sus1A, sus2A, sus2A]
                    r[1].K = 1000000

            sim = Simulation('Wmdirect', mdl, geom, self.rng)
            resSaver = self.__class__._addSaverToSim(sim, [Comp1])

            for i in range(self.nbRuns):
                sim.newRun()
                sim.cyt.SA.Count = 400
                sim.cyt.SB.Count = 400
                sim.cyt.Comp1[(sus1A,)*2 + (sus2C,)*2].Count = 500
                sim.run(10)

            allData.append(np.array(resSaver.data[:]))

        self._testData(allData)

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
    all_tests.append(unittest.makeSuite(StepsComplexReaction, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
