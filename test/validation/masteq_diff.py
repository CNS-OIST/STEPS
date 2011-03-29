import math
import numpy
import time 
from pylab import *

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as ssolv
import steps.utilities.meshio as meshio

from tol_funcs import *

### NOW   A+B-> B,  0->A (see Erban and Chapman, 2009)

########################################################################
SCALE = 1.0

KCST_f = 100e6*SCALE			# The reaction constant, degradation
KCST_b = (20.0e-10*SCALE) 		# The reaction constant, production

DCST_A = 20e-12
DCST_B = 20e-12

B0 = 1   # The number of B moleucles

DT = 0.1			# Sampling time-step
INT = 50000.1 		# Sim endtime

filename = 'cube_1_1_1_73tets.inp'

# A tolerance of 7.5% will fail <1% of the time
tolerance = 7.5/100

########################################################################

mdl  = smod.Model()

A = smod.Spec('A', mdl)
B = smod.Spec('B', mdl)

volsys = smod.Volsys('vsys',mdl)

diffA = smod.Diff('diffA', volsys, A, dcst = DCST_A)
diffB = smod.Diff('diffB', volsys, B, dcst = DCST_B)

# Production
R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [B], kcst = KCST_f)
R2 = smod.Reac('R2', volsys, lhs = [], rhs = [A], kcst = KCST_b)

geom = meshio.loadMesh('validation/meshes/'+filename)[0]

comp1 = sgeom.TmComp('comp1', geom, range(geom.ntets))
comp1.addVolsys('vsys')

rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295))

sim = ssolv.Tetexact(mdl, geom, rng)
sim.reset()

tpnts = numpy.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

res = numpy.zeros([ntpnts])
res_std1 = numpy.zeros([ntpnts])
res_std2 = numpy.zeros([ntpnts])

sim.reset()
sim.setCompCount('comp1', 'A', 0)
sim.setCompCount('comp1', 'B', B0)

b_time = time.time()
for t in range(0, ntpnts):
    sim.run(tpnts[t])
    res[t] = sim.getCompCount('comp1', 'A')


def fact(x): return (1 if x==0 else x * fact(x-1))

# Do cumulative count, but not comparing them all. 
# Don't get over 50 (I hope)
steps_n_res = numpy.zeros(50)
for r in res: steps_n_res[r]+=1
for s in range(50): steps_n_res[s] = steps_n_res[s]/ntpnts

passed= True
max_err = 0.0

k1 = KCST_f/6.022e23
k2 = KCST_b*6.022e23
v = comp1.getVol()*1.0e3 # litres

for m in range(5, 11):
    analy = (1.0/fact(m))*math.pow((k2*v*v)/(B0*k1), m)*math.exp(-((k2*v*v)/(k1*B0)))
    if not tolerable(steps_n_res[m], analy, tolerance):
        passed = False
    #if (abs(2*(steps_n_res[m]-analy)/(steps_n_res[m]+analy)) > max_err): max_err = abs(2*(steps_n_res[m]-analy)/(steps_n_res[m]+analy))
#print "Max error: ", max_err*100, "%"


