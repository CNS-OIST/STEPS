# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# This file is part of STEPS.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""
"""


from steps.geom import Geom
from steps.geom import Comp
from steps.model import Model
from steps.model import Reac
from steps.model import Spec
from steps.model import Volsys

import steps.rng
import steps.tetexact as stetexact


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def model():
    mdl = Model()
    x = Spec('X', mdl)
    y = Spec('Y', mdl)
    vsys = Volsys('main', mdl)
    rf = Reac('Rf', vsys, lhs=[x], rhs=[y])
    rb = Reac('Rb', vsys, lhs=[y], rhs=[x])
    rf.kcst = 1
    rb.kcst = 1
    return mdl


def geom():
    gm = Geom()
    c = Comp('comp', gm)
    c.volsys = set(['main'])
    c.vol = 1.6667e-21
    return gm
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class tetexact(object):


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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
