import math
import numpy
import time 
from pylab import *

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv

from tol_funcs import *

########################################################################

KCST_f = 100/(6.022e23*1.0e-15) 			# The reaction constant, production
KCST_b = 10			# The reaction constant, degradation
VOL = 1.0e-18

DT = 0.1			# Sampling time-step
INT = 200000.1 		# Sim endtime

# Tolerance for the comparison:
# In tests with good code <1% fail with tolerance of 1.5%
tolerance = 1.5/100  

########################################################################

mdl  = smod.Model()

A = smod.Spec('A', mdl)

volsys = smod.Volsys('vsys',mdl)

# Production
R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A], kcst = KCST_f)
R2 = smod.Reac('R2', volsys, lhs = [A], rhs = [], kcst = KCST_b)

geom = sgeom.Geom()

comp1 = sgeom.Comp('comp1', geom, VOL)
comp1.addVolsys('vsys')

rng = srng.create('mt19937', 1000)
rng.initialize(int(time.time()%4294967295))

sim = ssolv.Wmdirect(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res = numpy.zeros([ntpnts])

sim.reset()
sim.setCompCount('comp1', 'A', 0)

for t in xrange(0, ntpnts):
    sim.run(tpnts[t])
    res[t] = sim.getCompCount('comp1', 'A')

def fact(x): return (1 if x==0 else x * fact(x-1))

# Do cumulative count, but not comparing them all. 
# Don't get over 50 (I hope)
steps_n_res = numpy.zeros(50)
for r in res: steps_n_res[r]+=1
for s in range(50): steps_n_res[s] = steps_n_res[s]/ntpnts

passed = True
max_err = 0.0

k1 = KCST_b
k2 = KCST_f*(6.022e23*1.0e-15)

# Compare 5 to 15
for m in range(5, 16):
    analy = (1.0/fact(m))*math.pow((k2/k1), m)*math.exp(-(k2/k1))
    if not tolerable(steps_n_res[m], analy, tolerance):
        passed = False
    if (abs(2*(steps_n_res[m]-analy)/(steps_n_res[m]+analy)) > max_err): max_err = abs(2*(steps_n_res[m]-analy)/(steps_n_res[m]+analy))
    #print abs(2*(steps_n_res[m]-analy)/(steps_n_res[m]+analy))         
#print "Max error: ", max_err*100, "%"

