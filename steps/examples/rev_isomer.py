
from steps.model import Model
from steps.model import Reaction
from steps.model import Species
from steps.model import Volsys

from steps.geom import Geom
from steps.geom import Comp

import steps.rng

import steps.sim.wmdirect

def model():
    mdl = Model()
    x = Species('X', mdl)
    y = Species('Y', mdl)
    vsys = Volsys('main', mdl)
    rf = Reaction('Rf', vsys, lhs=[x], rhs=[y])
    rb = Reaction('Rb', vsys, lhs=[y], rhs=[x])
    rf.kf = 1
    rb.kf = 1
    return mdl


def geom():
    gm = Geom()
    c = Comp('comp', gm)
    c.volsys = set(['main'])
    c.volume = 1.6667e-21
    return gm
    

class wmdirect(object):

    def __init__(self):
        self.rng = steps.rng.create('mt19937', 512)
        self.model = model()
        self.geom = geom()
        self.sim = steps.sim.wmdirect.WMDirect(self.model, self.geom, self.rng)
    
    def iterate(self, seed, dt, maxt):
        self.rng.initialize(seed)
        self.sim.reset()
        self.sim.run(0)
        
        self.sim.setCompCount('comp', 'X', 500)
        self.sim.setCompCount('comp', 'Y', 0)
        self.sim.setCompVol('comp', 1.6667e-21)
        self.sim.setCompReacK('comp', 'Rf', 1)
        self.sim.setCompReacK('comp', 'Rb', 1)
        
        t = self.sim.time
        self.inter()
        times = []
        xs = []
        ys = []
        times.append(t)
        xs.append(self.sim.getCompCount('comp', 'X'))
        ys.append(self.sim.getCompCount('comp', 'Y'))
        while t < maxt:
            t = t + dt
            if t > maxt: t = maxt
            self.sim.run(t)
            self.inter()
            times.append(t)
            xs.append(self.sim.getCompCount('comp', 'X'))
            ys.append(self.sim.getCompCount('comp', 'Y'))
            
        return (times, xs, ys)
    
    def inter(self):
        pass
