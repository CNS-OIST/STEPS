####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

from steps import stepslib

# -----------------------------------------------------
# Steps classes are named here, and possibly extended
# -----------------------------------------------------
class Model(stepslib._py_Model)     : pass
class Volsys(stepslib._py_Volsys)   : pass
class Surfsys(stepslib._py_Surfsys) : pass
class Spec(stepslib._py_Spec)       : pass
class Reac(stepslib._py_Reac)       : pass
class SReac(stepslib._py_SReac)     : pass
class Diff(stepslib._py_Diff)       : pass
class Chan(stepslib._py_Chan)       : pass
class ChanState(stepslib._py_ChanState): pass
class OhmicCurr(stepslib._py_OhmicCurr): pass
class GHKcurr(stepslib._py_GHKcurr) : pass

class _py_VDepTrans(stepslib._py_VDepTrans):
    def __init__(self, id, surfsys, src_chan, dst_chan, **kwargs ):
        rate_f = kwargs['rate']
        vtable = _VoltageTable(rate_f, kwargs.get('vrange'))
        # Construct
        super(self.__class__, self).__init__(id, surfsys, src_chan, dst_chan, **vtable.as_dict())

class VDepSReac(stepslib._py_VDepSReac):
    def __init__(self, id, surfsys, **kwargs):
        # Get optional arguments that must be passed without change
        all_args = { hs:val for hs, val in kwargs.items() if hs in {'olhs', 'ilhs', 'slhs', 'irhs', 'srhs', 'orhs'} }
        # Process vrange
        rate_f = kwargs['k']
        vtable = _VoltageTable(rate_f, kwargs.get('vrange'))
        # Construct
        all_args.update(vtable.as_dict())
        super(self.__class__, self).__init__(id, surfsys, **all_args)


# -----------------------------------------------------
# Aux: Voltage table creator
# -----------------------------------------------------
class _VoltageTable:
    # These parameters are the hard-coded default voltage bounds
    vmin = -150.0e-3
    vmax = 100.0e-3
    dv = 1.0e-4

    def __init__(self, k_func, vrange=None):
        self.k_func = k_func
        if vrange:
            assert len(vrange)==3
            self.vmin = vrange[0]
            self.vmax = vrange[1]
            self.dv   = vrange[2]

        self.tablesize = int((self.vmax - self.vmin) / self.dv) + 1
        klist = [0.0] * self.tablesize

        #Calc table
        v = self.vmin
        for i in range(self.tablesize):
            klist[i] = k_func(v)
            v += self.dv

        #save klist as ktab
        self.ktab = klist

    #----------------
    def as_dict(self):
        """Return all important properties as a dictionary"""
        indexes = ['vmin', 'vmax', 'dv', 'ktab', 'tablesize']
        return { idx: getattr(self, idx)
                   for idx in indexes }
