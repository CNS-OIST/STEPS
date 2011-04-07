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

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import numpy
from pylab import *
import time

from tol_funcs import *

########################################################################

KCST = 5			# The reaction constant
N = 50				# Can set count or conc
VOL = 1.0e-18

NITER = 100000		# The number of iterations
DT = 0.1			# Sampling time-step
INT = 1.1			# Sim endtime

# Tolerance for the comparison:
# In test runs, with good code, < 1%  will fail with a 1.5% tolerance
tolerance = 2.0/100  

########################################################################

mdl  = smod.Model()

A = smod.Spec('A', mdl)
volsys = smod.Volsys('vsys',mdl)
R1 = smod.Reac('R1', volsys, lhs = [A], rhs = [], kcst = KCST)


geom = sgeom.Geom()
comp1 = sgeom.Comp('comp1', geom, VOL)
comp1.addVolsys('vsys')

rng = srng.create('mt19937', 1000)
rng.initialize(int(time.time()%4294967295))


sim = ssolv.Wmdirect(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res_m = numpy.zeros([NITER, ntpnts, 1])
res_std1 = numpy.zeros([ntpnts, 1])
res_std2 = numpy.zeros([ntpnts, 1])

for i in range (0, NITER):
	sim.reset()
	sim.setCompCount('comp1', 'A', N)

	for t in range(0, ntpnts):
		sim.run(tpnts[t])
		res_m[i, t, 0] = sim.getCompCount('comp1', 'A')

mean_res = numpy.mean(res_m, 0)
std_res = numpy.std(res_m, 0)

m_tol = 0
s_tol=0

passed = True
for i in range(ntpnts):
    if i == 0: continue
    analy = N*math.exp(-KCST*tpnts[i])
    std = math.pow((N*(math.exp(-KCST*tpnts[i]))*(1-(math.exp(-KCST*tpnts[i])))), 0.5)
    if not tolerable(analy, mean_res[i], tolerance):
        passed = False
    if not tolerable(std, std_res[i], tolerance): 
        passed = False


