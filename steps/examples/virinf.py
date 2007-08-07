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


"""A simple stochastic model for a viral infection.

Example usage: 
    >>> import steps.examples.virinf as vi
    >>> s = vi.wmdirect()
    >>> (times, gs, rs, vs, ps) = s.iterate(8482, 200.0/500.0, 200.0)
    >>> from pylab import *
    >>> plot(times, rs)
    >>> show()

Author: Zhi Xie <xiezhi@gmail.com>
"""


from steps.geom import Geom
from steps.geom import Comp
from steps.model import Model
from steps.model import Reac
from steps.model import Spec
from steps.model import Volsys

import steps.rng
import steps.sim.wmdirect

import numpy


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def model():
    mdl = Model()
    
    # Declare the species in the model.
    g = Spec('G', mdl)          # (G)enome
    r = Spec('R', mdl)          # m(R)NA
    v = Spec('V', mdl)          # (P)rotein
    p = Spec('P', mdl)          # (V)irus
    
    # Declare the reaction channels.
    vsys = Volsys('main', mdl)
    # R1: G ---> R
    r1 = Reac('R1', vsys, lhs=[g], rhs=[r])
    r1.kcst = 0.025
    # R2: R ---> *
    r2 = Reac('R2', vsys, lhs=[r], rhs=[])
    r2.kcst = 0.25
    # R3: R ---> G + R
    r3 = Reac('R3', vsys, lhs=[r], rhs=[g,r])
    r3.kcst = 1
    # R4: G + P ---> V
    r4 = Reac('R4', vsys, lhs=[g,p], rhs=[v])
    r4.kcst = 7.5e-6
    # R5: R ---> R + P
    r5 = Reac('R5', vsys, lhs=[r], rhs=[r,p])
    r5.kcst = 1000
    # R6: P ---> *
    r6 = Reac('R6', vsys, lhs=[p], rhs=[])
    r6.kcst = 1.99
    
    return mdl


def geom():
    gm = Geom()
    c = Comp('cell', gm)
    c.volsys = set(['main'])
    c.vol = 1.66666667e-27
    return gm
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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
        
        self.sim.setCompCount('cell', 'G', 0)
        self.sim.setCompCount('cell', 'R', 5)
        self.sim.setCompCount('cell', 'V', 0)
        self.sim.setCompCount('cell', 'P', 0)
        
        t = self.sim.time
        self.inter()
        nstops = int(round((maxt / dt) - 0.5)) + 2
        times = numpy.zeros((nstops,1))
        count = numpy.zeros((nstops,4))
        props = numpy.zeros((nstops,6))
        ix = 0
        times[ix] = 0.0
        count[ix,0] = self.sim.getCompCount('cell', 'G')
        count[ix,1] = self.sim.getCompCount('cell', 'R')
        count[ix,2] = self.sim.getCompCount('cell', 'V')
        count[ix,3] = self.sim.getCompCount('cell', 'P')
        props[ix,0] = self.sim.getCompReacA('cell', 'R1')
        props[ix,1] = self.sim.getCompReacA('cell', 'R2')
        props[ix,2] = self.sim.getCompReacA('cell', 'R3')
        props[ix,3] = self.sim.getCompReacA('cell', 'R4')
        props[ix,4] = self.sim.getCompReacA('cell', 'R5')
        props[ix,5] = self.sim.getCompReacA('cell', 'R6')
        while t < maxt:
            t = t + dt
            ix = ix + 1
            if t > maxt: t = maxt
            self.sim.run(t)
            self.inter()
            
            #times.append(t)
            #gs.append(self.sim.getCompCount('cell', 'G'))
            #rs.append(self.sim.getCompCount('cell', 'R'))
            #vs.append(self.sim.getCompCount('cell', 'V'))
            #ps.append(self.sim.getCompCount('cell', 'P'))
            times[ix] = t
            count[ix,0] = self.sim.getCompCount('cell', 'G')
            count[ix,1] = self.sim.getCompCount('cell', 'R')
            count[ix,2] = self.sim.getCompCount('cell', 'V')
            count[ix,3] = self.sim.getCompCount('cell', 'P')
            props[ix,0] = self.sim.getCompReacA('cell', 'R1')
            props[ix,1] = self.sim.getCompReacA('cell', 'R2')
            props[ix,2] = self.sim.getCompReacA('cell', 'R3')
            props[ix,3] = self.sim.getCompReacA('cell', 'R4')
            props[ix,4] = self.sim.getCompReacA('cell', 'R5')
            props[ix,5] = self.sim.getCompReacA('cell', 'R6')
        
        return (times, count, props)

    
    def inter(self):
        pass
        

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
