import steps.model as smodel
import steps.geom.core as sgeom
import steps.rng as srng
import steps.wmdirect as swmdirect

import numpy
from pylab import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

mdl = smodel.Model()

molA = smodel.Spec('molA', mdl)
molB = smodel.Spec('molB', mdl)
molC = smodel.Spec('molC', mdl)

volsys = smodel.Volsys('vsys', mdl)

kreac_f = smodel.Reac('kreac_f', volsys, lhs=['molA',molB], rhs=[molC])
kreac_f.kcst = 0.3e6
kreac_b = smodel.Reac('kreac_b', volsys, lhs=[molC], rhs=[molA,molB])
kreac_b.kcst = 0.7

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

mesh = sgeom.Container()

comp = sgeom.Comp('comp', mesh)
comp.addVolsys('vsys')
comp.vol = 1.6667e-21

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

r = srng.create('mt19937', 512)

sim = swmdirect.Solver(mdl, mesh, r)

r.initialize(23412)

res = numpy.zeros([100,2001,3])
tpnt = numpy.arange(0.0, 2.001, 0.001)

for i in xrange(0,100):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)
    
    for t in xrange(0,2001):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

res2 = numpy.mean(res, 0)

plot(tpnt, res2[:,0])
plot(tpnt, res2[:,1])
plot(tpnt, res2[:,2])
xlabel('Time (sec)')
ylabel('#molecules')
legend(('A','B','C'))
savefig('secondorderreaction02.svg')
show()

