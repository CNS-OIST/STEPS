
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

#########################
# Model setup
#########################

mdl = smodel.Model()

# Calcium
Ca = smodel.Spec('Ca', mdl)

# IP3
IP3 = smodel.Spec('IP3', mdl)

############### receptor state objects ###############
# receptor state: 'naive' state (no bound ligands)
R = smodel.Spec('R', mdl)

# receptor state: bound IP3
RIP3 = smodel.Spec('RIP3', mdl)

# receptor state: bound IP3 and Ca (open)
Ropen = smodel.Spec('Ropen', mdl)

# receptor state: Ca bound to one inactivation site
RCa = smodel.Spec('RCa', mdl)

# receptor state: Ca bound to two inactivation sites
R2Ca = smodel.Spec('R2Ca', mdl)

# receptor state: Ca bound to three inactivation sites
R3Ca = smodel.Spec('R3Ca', mdl)

# receptor state: Ca bound to four inactivation sites
R4Ca = smodel.Spec('R4Ca', mdl)
#######################################################

surfsys = smodel.Surfsys('ssys', mdl)

# The 'forward' binding reactions:

R_bind_IP3_f = smodel.SReac('R_bind_IP3_f', surfsys, olhs=[IP3], slhs=[R], srhs=[RIP3])

RIP3_bind_Ca_f = smodel.SReac('RIP3_bind_Ca_f', surfsys, olhs=[Ca], slhs=[RIP3], srhs = [Ropen])

R_bind_Ca_f = smodel.SReac('R_bind_Ca_f', surfsys, olhs=[Ca], slhs=[R], srhs=[RCa])

RCa_bind_Ca_f = smodel.SReac('RCa_bind_Ca_f', surfsys, olhs=[Ca], slhs=[RCa],srhs = [R2Ca])

R2Ca_bind_Ca_f = smodel.SReac('R2Ca_bind_Ca_f', surfsys, olhs=[Ca], slhs= [R2Ca], srhs = [R3Ca])

R3Ca_bind_Ca_f = smodel.SReac('R3Ca_bind_ca_f', surfsys, olhs=[Ca], slhs=[R3Ca], srhs=[R4Ca])

# The 'backward' unbinding reactions:

R_bind_IP3_b = smodel.SReac('R_bind_IP3_b', surfsys, slhs=[RIP3], orhs=[IP3], srhs=[R])

RIP3_bind_Ca_b = smodel.SReac('RIP3_bind_Ca_b', surfsys, slhs=[Ropen], orhs=[Ca], srhs=[RIP3])

R_bind_Ca_b = smodel.SReac('R_bind_Ca_b', surfsys, slhs=[RCa], orhs=[Ca], srhs=[R])

RCa_bind_Ca_b = smodel.SReac('RCa_bind_Ca_b', surfsys, slhs=[R2Ca], orhs=[Ca], srhs=[RCa])

R2Ca_bind_Ca_b = smodel.SReac('R2Ca_bind_Ca_b', surfsys, slhs=[R3Ca], orhs=[Ca], srhs= [R2Ca])

R3Ca_bind_Ca_b = smodel.SReac('R3Ca_bind_ca_b', surfsys, slhs=[R4Ca], orhs=[Ca], srhs=[R3Ca])

# Ca ions passing through open IP3R channel
R_Ca_channel_f = smodel.SReac('R_Ca_channel_f', surfsys, ilhs=[Ca], slhs=[Ropen], orhs=[Ca], srhs=[Ropen])

R_bind_IP3_f.setKcst(1000e6)
R_bind_IP3_b.setKcst(25800)
RIP3_bind_Ca_f.setKcst(8000e6)
RIP3_bind_Ca_b.setKcst(2000)
R_bind_Ca_f.setKcst(8.889e6)
R_bind_Ca_b.setKcst(5)
RCa_bind_Ca_f.setKcst(20e6)
RCa_bind_Ca_b.setKcst(10)
R2Ca_bind_Ca_f.setKcst(40e6)
R2Ca_bind_Ca_b.setKcst(15)
R3Ca_bind_Ca_f.setKcst(60e6)
R3Ca_bind_Ca_b.setKcst(20)
R_Ca_channel_f.setKcst(2e8)

#########################
# Geom setup
#########################

import steps.geom as swm
wmgeom = swm.Geom()

# Create the cytosol compartment
cyt = swm.Comp('cyt', wmgeom)
cyt.setVol(1.6572e-19)

# Create the Endoplasmic Reticulum compartment
ER = swm.Comp('ER', wmgeom, vol = 1.968e-20)

# ER is the 'inner' compartment, cyt is the 'outer' compartment
memb = swm.Patch('memb', wmgeom, ER, cyt)
memb.addSurfsys('ssys')
memb.setArea(0.4143e-12)

print 'Inner compartment to memb is', memb.getIComp().getID()

print 'Outer compartment to patch is', memb.getOComp().getID()

#########################
# RNG setup
#########################
import steps.rng as srng
r = srng.create('mt19937', 512)
r.initialize(7233)

#########################
# Solver setup
#########################
import steps.solver as ssolver
sim = ssolver.Wmdirect(mdl, wmgeom, r)

sim.setCompConc('cyt', 'Ca', 3.30657e-8)
sim.setCompCount('cyt', 'IP3', 6)
sim.setCompConc('ER', 'Ca', 150e-6)
sim.setCompClamped('ER', 'Ca', True)
sim.setPatchCount('memb', 'R', 160)

#########################
# Run simulation
#########################

sim.run(1e-3, 1e-4)