import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

import math
import time 
import numpy
from pylab import *

from tol_funcs import *

########################################################################

KCST = 5.0e6			# The reaction constant

CONCA = 1.0e-6
n = 2
CONCB = CONCA/n
VOL = 9.0e-18			#4.1533524410151228e-18

NITER = 1000			# The number of iterations
DT = 0.1			# Sampling time-step
INT = 1.1			# Sim endtime

# In test runs, with good code, <0.1% will fail with a tolerance of 1% 
tolerance = 1.0/100

########################################################################

mdl  = smod.Model()

A = smod.Spec('A', mdl)
B = smod.Spec('B', mdl)
C = smod.Spec('C', mdl)

volsys = smod.Volsys('vsys',mdl)

R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [C], kcst = KCST)

geom = sgeom.Geom()

comp1 = sgeom.Comp('comp1', geom, VOL)
comp1.addVolsys('vsys')

rng = srng.create('mt19937', 512)
import random
#rng.initialize(int(time.time()%4294967295))
import random
rng.initialize(int(random.random()*1000))

sim = ssolv.Wmdirect(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res_m = numpy.zeros([NITER, ntpnts, 3])


for i in range (0, NITER):
	sim.reset()
	sim.setCompConc('comp1', 'A', CONCA)
	sim.setCompConc('comp1', 'B', CONCB)

	for t in xrange(0, ntpnts):
		sim.run(tpnts[t])
		res_m[i, t, 0] = sim.getCompConc('comp1', 'A')
		res_m[i, t, 1] = sim.getCompConc('comp1', 'B')

mean_res = numpy.mean(res_m, 0)

lnBA = numpy.zeros(ntpnts)
lineAB = numpy.zeros(ntpnts)
C = CONCA-CONCB
passed  =True
max_err = 0.0
for i in range(ntpnts):
    A = mean_res[i][0]
    B = mean_res[i][1]
    lnBA[i] = math.log(B/A)
    lineAB[i] = math.log(CONCB/CONCA) -C*KCST*tpnts[i]

    if not tolerable(lnBA[i], lineAB[i], tolerance): passed = False
"""    if (abs(2*(lnBA[i]-lineAB[i])/(lnBA[i]+lineAB[i])) > max_err): max_err = abs(2*(lnBA[i]-lineAB[i])/(lnBA[i]+lineAB[i]))

print "max error", max_err*100.0, "%"
"""