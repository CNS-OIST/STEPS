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

res = numpy.zeros([100,12001,3])
tpnt = numpy.arange(0.0, 12.001, 0.001)

def run(i, tp1, tp2):
    for t in xrange(tp1,tp2):
        sim.run(tpnt[t])
        res[i,t,0] = sim.getCompCount('comp', 'molA')
        res[i,t,1] = sim.getCompCount('comp', 'molB')
        res[i,t,2] = sim.getCompCount('comp', 'molC')

for i in xrange(0,100):

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

res2 = numpy.mean(res, 0)

plot(tpnt, res2[:,0])
plot(tpnt, res2[:,1])
plot(tpnt, res2[:,2])
xlabel('Time (sec)')
ylabel('#molecules')
legend(('A','B','C'))
savefig('secondorderreaction08.svg')
show()
