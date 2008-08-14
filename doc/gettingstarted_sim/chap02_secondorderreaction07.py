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

res = numpy.zeros([12001,3])
tpnt = numpy.arange(0.0, 12.001, 0.001)

def run(tp1, tp2):
    for t in xrange(tp1,tp2):
        sim.run(tpnt[t])
        res[t,0] = sim.getCompCount('comp', 'molA')
        res[t,1] = sim.getCompCount('comp', 'molB')
        res[t,2] = sim.getCompCount('comp', 'molC')

sim.reset()
sim.setCompConc('comp', 'molA', 31.4e-6)
sim.setCompConc('comp', 'molB', 22.3e-6)

run(0,2001)
sim.setCompReacActive('comp', 'kreac_f', False)
run(2001,4001)
sim.setCompReacActive('comp', 'kreac_f', True)
run(4001,6001)
sim.setCompReacActive('comp', 'kreac_b', False)
run(6001,8001)
sim.setCompReacActive('comp', 'kreac_b', True)
run(8001,10001)
sim.setCompReacActive('comp', 'kreac_f', False)
sim.setCompReacActive('comp', 'kreac_b', False)
run(10001,12001)

plot(tpnt, res[:,0])
plot(tpnt, res[:,1])
plot(tpnt, res[:,2])
xlabel('Time (sec)')
ylabel('#molecules')
legend(('A','B','C'))
savefig('secondorderreaction07.svg')
show()
