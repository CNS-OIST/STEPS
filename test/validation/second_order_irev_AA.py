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

KCST = 50.0e6			# The reaction constant

CONCA = 20.0e-6
CONCB = CONCA

VOL = 9.0e-18			

NITER = 1000			# The number of iterations
DT = 0.1                # Sampling time-step
INT = 1.1               # Sim endtime

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
rng.initialize(int(time.time()%4294967295))

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

invA = numpy.zeros(ntpnts)
invB = numpy.zeros(ntpnts)
lineA  = numpy.zeros(ntpnts)
lineB = numpy.zeros(ntpnts)

max_err=0.0
passed = True
for i in range(ntpnts):
    invA[i] = (1.0/mean_res[i][0])
    invB[i] = (1.0/mean_res[i][1])
    lineA[i] = (1.0/CONCA +((tpnts[i]*KCST)))
    lineB[i] = (1.0/CONCB + ((tpnts[i]*KCST)))
    
    if not tolerable(invA[i], lineA[i], tolerance): passed = False
    if not tolerable(invB[i], lineB[i], tolerance): passed = False
"""
    if (abs(2*(invA[i]-lineA[i])/(invA[i]+lineA[i])) > max_err): max_err = abs(2*(invA[i]-lineA[i])/(invA[i]+lineA[i]))
    if (abs(2*(invB[i]-lineB[i])/(invB[i]+lineB[i])) > max_err): max_err = abs(2*(invB[i]-lineB[i])/(invB[i]+lineB[i]))
    
    
print "max error", max_err*100.0, "%"
"""