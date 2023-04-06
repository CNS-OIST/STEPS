####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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
class Model(stepslib._py_Model)     : 
    """
    Top-level container for the objects in a kinetic model.
    
    """
    pass


class Volsys(stepslib._py_Volsys)   : 
    """
    A container that groups reactions and diffusion rules involving Species located in a volume.
    """
    pass


class Surfsys(stepslib._py_Surfsys) : 
    """
    A container that groups reactions, diffusion rules and other transport mechanisms involving a Species located in a membrane.
    """
    pass


class Spec(stepslib._py_Spec)       : 
    """
    A chemical species which can be a reactant and/or product in reaction
    stoichiometry and/or associated with a diffusion rule, or other transport mechanisms.
    """
    pass


class Reac(stepslib._py_Reac)       : 
    """
    A kinetic reaction rule in a volume.
    
    The reaction rule is specified by:

    * Species on the left hand side of the reaction: the reactants (lhs).
    * Species on the right hand side of the reaction: the products (rhs).
    * Rate constant for the reaction, supplied in s.i. units (kcst).
    """
    pass


class SReac(stepslib._py_SReac)     : 
    """
    A reaction rule where at least one reactant is embedded in a surface.
    
    In a surface reaction, the species can be classified as:

    * Reactants on the left hand side of the reaction:
        * Species in the 'outer' compartment (olhs)
        * or Species in the 'inner' compartment (ilhs)
        * Species on the surface (slhs)
    * Reactants on right hand side of the reaction:
        * Species in the 'outer' compartment (orhs)
        * Species in the 'inner' compartment (irhs)
        * Species on the surface (srhs)
    * The reaction rate is defined by kcst, supplied in s.i. units.
    """
    pass


class Diff(stepslib._py_Diff)       : 
    """
    A diffusion rule for a chemical species in a volume or on a surface.
    
    A diffusion rule is described by:

    * Species to which the diffusion rule applies (lig).
    * Diffusion constant (dcst) specified in s.i. units.
    """
    pass


class Chan(stepslib._py_Chan)       : 
    """
    A channel object which groups a set of channel states. Whilst not involved
    directly in reactions, currents, transitions etc it provides necessary
    grouping for channel states.
    """
    pass


class ChanState(stepslib._py_ChanState): 
    """
    A channel state object which may be involved in channel events such as
    currents, voltage-dependent transitions and any other steps.model.Spec events,
    such as surface reactions and diffusion.
    """
    pass


class OhmicCurr(stepslib._py_OhmicCurr): 
    """
    An ohmic current object that describes an ohmic current through a channel in
    a particular conducting state.
    
    An ohmic current object is a simple current based on a single, fixed value
    for single-channel conductance and a constant reversal potential. An ohmic current
    does not result in movement of ions between compartments, but simply
    contributes a continous current to the EField solver for every conducting channel
    on a membrane surface that contains the ohmic current.
    """
    pass


class GHKcurr(stepslib._py_GHKcurr) : 
    """
    A current object based on the Goldman-Hodgkin-Katz flux equation,
    that describes a current through a channel in a particular permeable state.
    
    Each GHK current in the simulation is solved within the SSA with a rate determined
    from the simulation state, i.e. membrane potential, 'outer' and 'inner'
    concentration of the ion and temperature (which is fixed), and constant parameters
    ion valence and permeability.
    For a permeabilty to be found the user must supply
    information from a channel measurement:

    * the single-channel conductance
    * the potential
    * the temparature (this may be different from the STEPS simulation temperature)
    * the 'outer' concentration of the ion
    * the 'inner' concentration of the ion
    
    This information must be set with the member function setPInfo()
    
    A GHK current involves, optionally, a real transfer of ions between comparments
    separated by a membrane in the STEPS simulation. It is important that
    these ions implement diffusion if compartments are not well-mixed.
    """
    pass



class _py_VDepTrans(stepslib._py_VDepTrans):
    def __init__(self, id, surfsys, src_chan, dst_chan, **kwargs ):
        rate_f = kwargs['rate']
        vtable = _VoltageTable(rate_f, kwargs.get('vrange'))
        # Construct
        super(self.__class__, self).__init__(id, surfsys, src_chan, dst_chan, **vtable.as_dict())

class VDepSReac(stepslib._py_VDepSReac):
    """
    A voltage-dependent reaction involving any species or channel state.
    
    A voltage-dependent reaction is similar to a surface reaction, except the
    reaction parameter is voltage-dependent.
    
    A function object must be passed to the constructor returning the reaction
    'constant' as a function of voltage. This function is used to fill a table of
    reaction 'constants' at a range of voltage values (default: -150mV to 100mV,
    in steps of 0.1mV) These default values may be altered with optional argument
    'vrange'.
    """
    def __init__(self, id, surfsys, **kwargs):
        """
        Construction::
        
            vdepsreac = steps.model.VDepSReac(id, surfsys,
            ilhs = [ ], olhs = [ ], slhs = [ ],
            irhs = [ ], orhs = [ ], srhs = [ ],
            k = <function>,
            vrange = [-150.0e-3, 100.0e-3, 1.0e-4] )
        
        Construct a voltage-dependent reaction object with identifier string id
        and assign surfsys as the parent surface system. A list of
        left hand reactants are assigned with ilhs, olhs and slhs
        (default for each is an empty list). A list of right hand side
        products are assigned with irhs, orhs and srhs (default for each
        is an empty list).
        A function that returns the kinetic reaction 'constant' in ordinary, Molar
        units at any voltage (in volts) is supplied with argument k.
        A 'voltage range' over which to calculate the reaction rate is
        optionally provided by the argument vrange (in Volts) as a list:
        [minimum voltage, maximum voltage, voltage step]
        
        Arguments:
            * string id
            * steps.model.Surfsys surfsys
            * list(steps.model.Spec) ilhs (default = [ ])
            * list(steps.model.Spec) olhs (default = [ ])
            * list(steps.model.Spec) slhs (default = [ ])
            * list(steps.model.Spec) irhs (default = [ ])
            * list(steps.model.Spec) orhs (default = [ ])
            * list(steps.model.Spec) srhs (default = [ ])
            * function k
            * list vrange (default = [-150.0e-3, 100.0e-3, 1.0e-4])
        """
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
