# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2011 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Example: Well-mixed reaction system
# http://steps.sourceforge.net/manual/well_mixed.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model as smodel

mdl = smodel.Model()

molA = smodel.Spec('molA', mdl)
molB = smodel.Spec('molB', mdl)
molC = smodel.Spec('molC', mdl)

vsys = smodel.Volsys('vsys', mdl)

kreac_f = smodel.Reac('kreac_f', vsys, lhs=[molA, molB], rhs=[molC], kcst=0.3e6)
kreac_b = smodel.Reac('kreac_b', vsys, lhs=[molC], rhs=[molA, molB])
kreac_b.kcst = 0.7

import steps.geom as swm

wmgeom = swm.Geom()

comp = swm.Comp('comp', wmgeom)
comp.addVolsys('vsys')
comp.setVol(1.6667e-21)

import steps.rng as srng
r = srng.create('mt19937', 256)
r.initialize(23412)

import steps.solver as ssolver
sim = ssolver.Wmdirect(mdl, wmgeom, r)

sim.reset()

sim.setCompConc('comp', 'molA', 31.4e-6)
sim.setCompConc('comp', 'molB', 22.3e-6)

import numpy
tpnt = numpy.arange(0.0, 2.001, 0.001)
res = numpy.zeros([2001,3])

for t in range(0,2001):
    sim.run(tpnt[t])
    res[t,0] = sim.getCompCount('comp', 'molA')
    res[t,1] = sim.getCompCount('comp', 'molB')
    res[t,2] = sim.getCompCount('comp', 'molC')

import pylab

# Plot number of molecules of 'molA' over the time range:
pylab.plot(tpnt, res[:,0], label = 'A')
# Plot number of molecules of 'molB' over the time range:
pylab.plot(tpnt, res[:,1], label = 'B')
# Plot number of molecules of 'molC' over the time range:
pylab.plot(tpnt, res[:,2], label = 'C')

pylab.xlabel('Time (sec)')
pylab.ylabel('#molecules')
pylab.legend()
pylab.show()

def plotres(res):
    # Plot number of molecules of 'molA' over the time range:
    pylab.plot(tpnt, res[:,0], label = 'A')
    # Plot number of molecules of 'molB' over the time range:
    pylab.plot(tpnt, res[:,1], label = 'B')
    # Plot number of molecules of 'molC' over the time range:
    pylab.plot(tpnt, res[:,2], label = 'C')

    pylab.xlabel('Time (sec)')
    pylab.ylabel('#molecules')
    pylab.legend()
    pylab.show()    

NITER = 100
res = numpy.zeros([NITER,2001,3])
tpnt = numpy.arange(0.0, 2.001, 0.001)

for i in range(0,NITER):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)

    for t in range(0,2001):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

res_mean = numpy.mean(res, 0)

plotres(res_mean)

for i in range(NITER):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)

    for t in range(0,1001):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

    # Add 10 molecules of species A
    sim.setCompCount('comp', 'molA', sim.getCompCount('comp', 'molA') + 10)
    for t in range(1001,2001):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

res_mean = numpy.mean(res, 0)
plotres(res_mean)

NITER = 1
res = numpy.zeros([NITER,2001,3])

for i in range(NITER):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)

    for t in range(0,101):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

    sim.setCompClamped('comp', 'molA', True)

    for t in range(101,601):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

    sim.setCompClamped('comp', 'molA', False)

    for t in range(601,2001):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')
    
res_mean = numpy.mean(res, 0)
plotres(res_mean)

NITER=100

def run(i, tp1, tp2):
    for t in range(tp1,tp2):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

res = numpy.zeros([NITER,12001,3])
tpnt = numpy.arange(0.0, 12.001, 0.001)

for i in range(NITER):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)
    run(i,0,2001)
    sim.setCompReacActive('comp', 'kreac_f', False)
    run(i,2001,4001)
    sim.setCompReacActive('comp', 'kreac_f', True)
    run(i,4001,6001)
    sim.setCompReacActive('comp', 'kreac_b', False)
    run(i,6001,8001)
    sim.setCompReacActive('comp', 'kreac_b', True)
    run(i,8001,10001)
    sim.setCompReacActive('comp', 'kreac_f', False)
    sim.setCompReacActive('comp', 'kreac_b', False)
    run(i,10001,12001)

res_mean = numpy.mean(res, 0)
plotres(res_mean)

