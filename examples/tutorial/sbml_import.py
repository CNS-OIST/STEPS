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

# Example: SBML import
# http://steps.sourceforge.net/manual/sbml_importer.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.rng as srng
import steps.solver as ssolver
import steps.utilities.sbml as ssbml
import numpy
from pylab import *

def runSBMLmod(sbmlFile, time_end, time_dt, iter_n = 1, solver = 'Wmdirect', specs_plot = {}, vol_units = 1.0e-3, vol_def = False):
    iSbml = ssbml.Interface(sbmlFile, volunits_def = vol_units, volume_def = vol_def )
    mdl = iSbml.getModel()
    mesh = iSbml.getGeom()
    comp_specs = {}
    if not specs_plot: got_sp = False
    else: got_sp = True
    for comp in mesh.getAllComps():
        comp_specs[comp.getID()] = []
        volsys_strings = comp.getVolsys()
        for vsys_str in volsys_strings:
            vsys = mdl.getVolsys(vsys_str)
            specs = vsys.getAllSpecs()
            for spec in specs:
                comp_specs[comp.getID()].append(spec.getID())
                if (got_sp == False): specs_plot.update({spec.getID():''})
    comp_specs_n = 0
    for comp in comp_specs:
        comp_specs[comp].sort()
        comp_specs_n += len(comp_specs[comp])
    time_pnts = numpy.arange(0.0, time_end, time_dt)
    points_n = int(round(time_end/time_dt))
    if (len(time_pnts)  == (points_n + 1)): time_pnts = numpy.delete(time_pnts, len(time_pnts)-1)
    res = numpy.zeros([iter_n, points_n, comp_specs_n])
    r = srng.create('mt19937', 256)
    r.initialize(7233)
    if (solver == 'Wmdirect'):
        sim = ssolver.Wmdirect(mdl, mesh, r)
    elif (solver == 'Wmrk4'):
        sim = ssolver.Wmrk4(mdl, mesh, r)
        sim.setDT(0.0001)
    else:
        raise NameError("Unsupported solver. SBML importer accepts well-mixed solvers 'Wmrk4' and 'Wmdirect'")
    for it in range (0, iter_n):
        sim.reset()
        iSbml.reset()
        iSbml.setupSim(sim)
        for tp in range(0, points_n):
            sim.run(time_pnts[tp])
            i = 0
            for comp in comp_specs:
                for spec in comp_specs[comp]:
                    res[it, tp, i] = sim.getCompConc(comp, spec)
                    i+=1
            iSbml.updateSim(sim, time_dt)
    mean_res = numpy.mean(res, 0)
    i=0
    for comp in comp_specs:
        for spec in comp_specs[comp]:
            if (spec in specs_plot):
                if (specs_plot[spec]): plot(time_pnts, mean_res[:,i], label = spec+ ", "+ comp, color = specs_plot[spec])
                else: plot(time_pnts, mean_res[:,i], label = spec+ ", "+ comp)
            i+=1
    legend(loc = 'best', numpoints=1)
    xlabel('Time (s)')
    ylabel('Conc (M)')
    show()

print "Running Biomodel 98 in deterministic solver..."
runSBMLmod('biomodel/BIOMD0000000098.xml', 10, 0.001, solver='Wmrk4',specs_plot={'Y':'blue', 'Z':'red'})

print "\nRunning Biomodel 98 in well-mixed stochastic solver in reduced volume..."
runSBMLmod('biomodel/BIOMD0000000098.xml', 10, 0.001, solver='Wmdirect',specs_plot={'Y':'blue', 'Z':'red'}, vol_units = 1.0e-3, vol_def = 1.0e-18)

