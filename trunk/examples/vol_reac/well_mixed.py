
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
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

import steps.model as smodel
import steps.geom as swm
import steps.rng as srng
import steps.solver as ssolver

##############################
# Model Setup
##############################

mdl = smodel.Model()

molA = smodel.Spec('molA', mdl)
molB = smodel.Spec('molB', mdl)
molC = smodel.Spec('molC', mdl)

vsys = smodel.Volsys('vsys', mdl)

kreac_f = smodel.Reac('kreac_f', vsys, lhs=[molA,molB], rhs=[molC], kcst=0.3e6)
kreac_b = smodel.Reac('kreac_b', vsys, lhs=[molC], rhs=[molA,molB])
kreac_b.kcst = 0.7

##############################
# Geom Setup
##############################
wmgeom = swm.Geom()

comp = swm.Comp('comp', wmgeom)
comp.addVolsys('vsys')
comp.setVol(1.6667e-21)

##############################
# RNG Setup
##############################
r = srng.create('mt19937', 256)
r.initialize(23412)

##############################
# Well Mixed DirectSSA
##############################
sim = ssolver.Wmdirect(mdl, wmgeom, r)

sim.reset()

sim.setCompConc('comp', 'molA', 31.4e-6)
sim.setCompConc('comp', 'molB', 22.3e-6)

###################################
# Checkpoint manually
###################################
print "Manually checkpoint"
sim.checkpoint("manual.checkpoint")
print "Sim Time: %e" % sim.getTime()
print "Sim Steps: %e" % sim.getNSteps()

###################################
# Run Simulation with checkpointing
###################################
print "Run to 2e-3s, checkpoint for every 1e-4s"
sim.run(2e-3, 1e-4)
print "Sim Time: %e" % sim.getTime()
print "Sim Steps: %e" % sim.getNSteps()

# or advance
print "advance 1e-3s checkpoint for every 1e-4s"
sim.advance(1e-3, 5e-4)
print "Sim Time: %e" % sim.getTime()
print "Sim Steps: %e" % sim.getNSteps()

# or step, without checkpointing
print "Run 1 SSA step"
sim.step()
print "Sim Time: %e" % sim.getTime()
print "Sim Steps: %e" % sim.getNSteps()

###################################
# Restore from checkpoint file
###################################
print "Restore from file"
sim.restore("manual.checkpoint")
print "Sim Time: %e" % sim.getTime()
print "Sim Steps: %e" % sim.getNSteps()

###################################
# RK4 solver
###################################
print "Well Mixed RK4"
rk4sim = ssolver.Wmrk4(mdl, wmgeom, r)

rk4sim.reset()
rk4sim.setDT(1e-4)
rk4sim.setCompConc('comp', 'molA', 31.4e-6)
rk4sim.setCompConc('comp', 'molB', 22.3e-6)

# Notice: No checkpointing for RK4 yet
rk4sim.run(1e-3)
print "Sim Time: %e" % rk4sim.getTime()