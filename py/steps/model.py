# -*- coding: utf-8 -*-

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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This file is the user-interface file for all model objects in steps.
# All objects are directly derived from the corresponding swig objects.
# Model container object is owned by Python
# All other objects are owned by c++ and Model container is responsible for 
# all the cleaning-up of these objects (see cpp/model/model.cpp class 
# destructor).
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

try:
    from . import steps_swig_numpy as steps_swig
    import _steps_swig_numpy as _steps_swig
except:
    from . import steps_swig
    import _steps_swig

import steps

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Model(steps_swig.Model) : 
    """
    Top-level container for the objects in a kinetic model.
    
    """
    
    def __init__(self, *args):
        """
        Construction::
        
            m = steps.model.Model()
            
        Create a model container object.
            
        Arguments: 
            None
        """
        
        this = _steps_swig.new_Model(*args)
        try: self.this.append(this)
        except: self.this = this
        # let the model object do the cleaning-up
        self.thisown = 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Spec(steps_swig.Spec) :

    """
    A chemical species which can be a reactant and/or product in reaction 
    stoichiometry and/or associated with a diffusion rule, or other transport mechanisms. 
    """
    def __init__(self, *args): 
        """        
        Construction::
        
            s = steps.model.Spec(id, mdl)
            
        Create a species object with identifier string id and assign the 
        object mdl as its parent model.
            
        Arguments: 
            * string id
            * steps.model.Model mdl
        """
        this = _steps_swig.new_Spec(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Spec_setID
        self.__swig_getmethods__["id"] = _steps_swig.Spec_getID
        self.__swig_getmethods__["model"] = _steps_swig.Spec_getModel
    id = steps_swig._swig_property(_steps_swig.Spec_getID, _steps_swig.Spec_setID)
    """Identifier string of the species."""
    model = steps_swig._swig_property(_steps_swig.Spec_getModel)
    """Reference to parent model."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Chan(steps_swig.Chan):
    
    """
        A channel object which groups a set of channel states. Whilst not involved 
        directly in reactions, currents, transitions etc it provides necessary 
        grouping for channel states.
        """
    def __init__(self, *args): 
        """
            Construction:
            
            chan = steps.model.Chan(id, mdl)
            
            Create a channel object with identifier string id and assign the object
            mdl as its parent model.
            
            Arguments: 
            * string id
            * steps.model.Model mdl        
            """
        this = _steps_swig.new_Chan(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Chan_setID
        self.__swig_getmethods__["id"] = _steps_swig.Chan_getID
        self.__swig_getmethods__["model"] = _steps_swig.Chan_getModel
    id = steps_swig._swig_property(_steps_swig.Chan_getID, _steps_swig.Chan_setID)
    """ Identifier string of the channel. """
    model = steps_swig._swig_property(_steps_swig.Chan_getModel)  
    """ Refernece to parent model. """

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class ChanState(steps_swig.ChanState):
    
    """
        A channel state object which may be involved in channel events such as 
        currents, voltage-dependent transitions and any other steps.model.Spec events,
        such as surface reactions and diffusion.
        """
    def __init__(self, *args): 
        """
            Construction:
            
            chanstate = steps.model.ChanState(id, mdl, chan)
            
            Create a channel state object with identifier string id, assign the object
            mdl as its parent model. This object describes one of the possible states
            of channel object chan.
            
            Arguments: 
            * string id
            * steps.model.Model mdl       
            * steps.model.Chan chan
            """
        this = _steps_swig.new_ChanState(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_getmethods__["chan"] = _steps_swig.ChanState_getChan  
    chan = steps_swig._swig_property(_steps_swig.ChanState_getChan)
    """ Reference to parent channel. """

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Surfsys(steps_swig.Surfsys):
    """
    A container that groups reactions, diffusion rules and other transport mechanisms involving a Species located in a membrane.
    """
    def __init__(self, *args): 
        """
        Construction::
        
            s = steps.model.Surfsys(id, mdl)
            
        Construct a surface system object with identifier string id and assign 
        the object mdl as its parent model.
        
        Arguments: 
            * string id
            * steps.model.Model mdl
        """
        this = _steps_swig.new_Surfsys(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Surfsys_setID
        self.__swig_getmethods__["id"] = _steps_swig.Surfsys_getID
        self.__swig_getmethods__["model"] = _steps_swig.Surfsys_getModel
    id = steps_swig._swig_property(_steps_swig.Surfsys_getID, _steps_swig.Surfsys_setID)
    """Identifier string of the surface system."""
    model = steps_swig._swig_property(_steps_swig.Surfsys_getModel)
    """Reference to parent model."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Volsys(steps_swig.Volsys) :
    """
    A container that groups reactions and diffusion rules involving Species located in a volume.
    """
    def __init__(self, *args): 
        """
        Construction::
        
            v = steps.model.Volsys(id, mdl)
            
        Construct a volume system object with identifier string id and assign 
        the object mdl as its parent model.
        
        Arguments: 
            * string id
            * steps.model.Model mdl
        """
        this = _steps_swig.new_Volsys(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Volsys_setID
        self.__swig_getmethods__["id"] = _steps_swig.Volsys_getID
        self.__swig_getmethods__["model"] = _steps_swig.Volsys_getModel
    id = steps_swig._swig_property(_steps_swig.Volsys_getID, _steps_swig.Volsys_setID)
    """Identifier string of the volume system."""
    model = steps_swig._swig_property(_steps_swig.Volsys_getModel)
    """Reference to parent model."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Diff(steps_swig.Diff) :
    """
    A diffusion rule for a chemical species in a volume or on a surface.
    
    A diffusion rule is described by:
    * Species to which the diffusion rule applies (lig).
    * Diffusion constant (dcst) specified in s.i. units.
    """
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            diff = steps.model.Diff(id, volsys, lig, dcst = 0.0)
            
        Construct a diffusion rule object with identifier string id applied to species 
        lig and assign volsys as the parent volume system. Diffusion constant is 
        set by dcst.
        
        Arguments: 
            * string id
            * steps.model.Volsys volsys
            * steps.model.Spec lig
            * float dcst (default = 0.0)
        
        
        Construction::
            
            diff = steps.model.Diff(id, surfsys, lig, dcst = 0.0)
            
        Construct a diffusion rule object with identifier string id applied to species 
        lig and assign surfsys as the parent surface system. Diffusion constant is 
        set by dcst.
            
        Arguments: 
            * string id
            * steps.model.Surfsys surfsys
            * steps.model.Spec lig
            * float dcst (default = 0.0)
        """
        
        this = _steps_swig.new_Diff(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Diff_setID
        self.__swig_getmethods__["id"] = _steps_swig.Diff_getID
        self.__swig_getmethods__["model"] = _steps_swig.Diff_getModel
        self.__swig_getmethods__["volsys"] = _steps_swig.Diff_getVolsys
        self.__swig_getmethods__["surfsys"] = _steps_swig.Diff_getSurfsys
        self.__swig_setmethods__["lig"] = _steps_swig.Diff_setLig
        self.__swig_getmethods__["lig"] = _steps_swig.Diff_getLig
        self.__swig_setmethods__["dcst"] = _steps_swig.Diff_setDcst
        self.__swig_getmethods__["dcst"] = _steps_swig.Diff_getDcst
    id = steps_swig._swig_property(_steps_swig.Diff_getID, _steps_swig.Diff_setID)
    """Identifier string of the diffusion."""
    model = steps_swig._swig_property(_steps_swig.Diff_getModel)
    """Reference to parent model."""
    volsys = steps_swig._swig_property(_steps_swig.Diff_getVolsys)
    """Reference to parent volume system."""
    surfsys = steps_swig._swig_property(_steps_swig.Diff_getSurfsys)
    """Reference to parent surface system."""
    dcst = steps_swig._swig_property(_steps_swig.Diff_getDcst, _steps_swig.Diff_setDcst)
    """Diffusion constant."""
    lig = steps_swig._swig_property(_steps_swig.Diff_getLig, _steps_swig.Diff_setLig)
    """Reference to diffusion species."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Reac(steps_swig.Reac) :
    """
    A kinetic reaction rule in a volume.
    
    The reaction rule is specified by:
    * Species on the left hand side of the reaction: the reactants (lhs).
    * Species on the right hand side of the reaction: the products (rhs).
    * Rate constant for the reaction, supplied in s.i. units (kcst).
    """
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            reac = steps.model.Reac(id, volsys, lhs = [ ], rhs = [ ], kcst = 0.0)
            
        Construct a reaction rule object with identifier string id and assign 
        volsys as the parent volume system. A list of left hand side reactants 
        may be assigned with lhs, whilst a list of right hand side products may 
        be assigned with rhs, the kinetic reaction rate constant is set by kcst.
        
        Arguments: 
            * string id
            * steps.model.Volsys volsys
            * list(steps.model.Spec) lhs (default = [ ]) 
            * list(steps.model.Spec) rhs (default = [ ]) 
            * float kcst (default = 0.0)
        """
        this = _steps_swig.new_Reac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Reac_setID
        self.__swig_getmethods__["id"] = _steps_swig.Reac_getID
        self.__swig_getmethods__["model"] = _steps_swig.Reac_getModel
        self.__swig_getmethods__["volsys"] = _steps_swig.Reac_getVolsys
        self.__swig_setmethods__["kcst"] = _steps_swig.Reac_setKcst
        self.__swig_getmethods__["kcst"] = _steps_swig.Reac_getKcst
        self.__swig_setmethods__["lhs"] = _steps_swig.Reac_setLHS
        self.__swig_getmethods__["lhs"] = _steps_swig.Reac_getLHS
        self.__swig_setmethods__["rhs"] = _steps_swig.Reac_setRHS
        self.__swig_getmethods__["rhs"] = _steps_swig.Reac_getRHS
        self.__swig_getmethods__["order"] = _steps_swig.Reac_getOrder
    id = steps_swig._swig_property(_steps_swig.Reac_getID, _steps_swig.Reac_setID)
    """Identifier string of the reaction."""
    model = steps_swig._swig_property(_steps_swig.Reac_getModel)
    """Reference to parent model."""
    volsys = steps_swig._swig_property(_steps_swig.Reac_getVolsys)
    """Reference to parent volume system."""
    kcst = steps_swig._swig_property(_steps_swig.Reac_getKcst, _steps_swig.Reac_setKcst)
    """Reaction constant."""
    lhs = steps_swig._swig_property(_steps_swig.Reac_getLHS, _steps_swig.Reac_setLHS)
    """Left hand side reactants."""
    rhs = steps_swig._swig_property(_steps_swig.Reac_getRHS, _steps_swig.Reac_setRHS)
    """Right hand side reactants."""
    order = steps_swig._swig_property(_steps_swig.Reac_getOrder)
    """Order of the reaction."""
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class SReac(steps_swig.SReac) :
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
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            sreac = steps.model.SReac(id, surfsys, 
                                    ilhs = [ ], olhs = [ ], slhs = [ ],
                                    irhs = [ ], orhs = [ ], srhs = [ ],
                                    kcst = 0.0)
            
        Construct a surface reaction rule object with identifier string 
        id and assign surfsys as the parent surface system. A list of 
        left hand reactants are assigned with ilhs, olhs and slhs 
        (default for each is an empty list). A list of right hand side 
        products are assigned with irhs, orhs and srhs (default for each 
        is an empty list). The kinetic reaction rate constant is set with kcst.
        
        
        Arguments: 
            * string id
            * steps.model.Surfsys surfsys
            * list(steps.model.Spec) ilhs (default = [ ]) 
            * list(steps.model.Spec) olhs (default = [ ])
            * list(steps.model.Spec) slhs (default = [ ])  
            * list(steps.model.Spec) irhs (default = [ ])
            * list(steps.model.Spec) orhs (default = [ ])
            * list(steps.model.Spec) srhs (default = [ ])
            * float kcst (default = 0.0)
        """
        this = _steps_swig.new_SReac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.SReac_setID
        self.__swig_getmethods__["id"] = _steps_swig.SReac_getID
        self.__swig_getmethods__["model"] = _steps_swig.SReac_getModel
        self.__swig_getmethods__["surfsys"] = _steps_swig.SReac_getSurfsys
        self.__swig_setmethods__["kcst"] = _steps_swig.SReac_setKcst
        self.__swig_getmethods__["kcst"] = _steps_swig.SReac_getKcst
        self.__swig_getmethods__["order"] = _steps_swig.SReac_getOrder
        self.__swig_getmethods__["outer"] = _steps_swig.SReac_getOuter
        self.__swig_setmethods__["olhs"] = _steps_swig.SReac_setOLHS
        self.__swig_getmethods__["olhs"] = _steps_swig.SReac_getOLHS
        self.__swig_setmethods__["ilhs"] = _steps_swig.SReac_setILHS
        self.__swig_getmethods__["ilhs"] = _steps_swig.SReac_getILHS
        self.__swig_setmethods__["slhs"] = _steps_swig.SReac_setSLHS
        self.__swig_getmethods__["slhs"] = _steps_swig.SReac_getSLHS
        self.__swig_setmethods__["irhs"] = _steps_swig.SReac_setIRHS
        self.__swig_getmethods__["irhs"] = _steps_swig.SReac_getIRHS
        self.__swig_setmethods__["srhs"] = _steps_swig.SReac_setSRHS
        self.__swig_getmethods__["srhs"] = _steps_swig.SReac_getSRHS
        self.__swig_setmethods__["orhs"] = _steps_swig.SReac_setORHS
        self.__swig_getmethods__["orhs"] = _steps_swig.SReac_getORHS
    id = steps_swig._swig_property(_steps_swig.SReac_getID, _steps_swig.SReac_setID)
    """Identifier string of the surface reaction."""
    model = steps_swig._swig_property(_steps_swig.SReac_getModel)
    """Reference to parent model."""
    surfsys = steps_swig._swig_property(_steps_swig.SReac_getSurfsys)
    """Reference to parent surface system."""
    kcst = steps_swig._swig_property(_steps_swig.SReac_getKcst, _steps_swig.SReac_setKcst)
    """Reaction constant."""
    outer = steps_swig._swig_property(_steps_swig.SReac_getOuter)
    """Obsolete"""
    olhs = steps_swig._swig_property(_steps_swig.SReac_getOLHS, _steps_swig.SReac_setOLHS)
    """Left hand side reactants in outer compartment."""
    ilhs = steps_swig._swig_property(_steps_swig.SReac_getILHS, _steps_swig.SReac_setILHS)
    """Left hand side reactants in inner compartment."""
    slhs = steps_swig._swig_property(_steps_swig.SReac_getSLHS, _steps_swig.SReac_setSLHS)
    """Left hand side reactants on surface."""
    irhs = steps_swig._swig_property(_steps_swig.SReac_getIRHS, _steps_swig.SReac_setIRHS)
    """Right hand side reactants in inner compartment."""
    srhs = steps_swig._swig_property(_steps_swig.SReac_getSRHS, _steps_swig.SReac_setSRHS)
    """Right hand side reactants on surface."""
    orhs = steps_swig._swig_property(_steps_swig.SReac_getORHS, _steps_swig.SReac_setORHS)
    """Right hand side reactants in outer compartment."""
    order = steps_swig._swig_property(_steps_swig.SReac_getOrder)
    """Order of the reaction."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class VDepTrans(steps_swig.VDepTrans) :
    
    """
        A voltage-dependent transition involving channel states embedded in a surface. 
        
        A voltage-dependent transition is a first-order Markov transition from one
        channel state to another. Thus two channel states descibe a voltage-dependent
        transition:
        The 'source' channel. This is the left-hand side of the reaction- the reactant.
        The 'destination' channel. This is the righ-hand side of the reaction- the product.
        
        A function object must be passed to the constructor returning the transition
        rate as a function of voltage. This function is used to fill a table of 
        transition rates at hard-coded voltage values (currently -150mV to 100mV, 
        in steps of 0.1mV) (These hard-coded values may be altered in the surface system
        in the future).
        """
    
    def __init__(self, *args, **kwargs): 
        """
            Construction::
            
            vdeptrans = steps.model.VDepTrans(id, surfsys, src, dst, rate = <function>) 
            
            Construct a voltage-dependent transition object with identifier string id
            and assign surfsys as the parent surface system. The 'source' 
            channel state is assigned with src and the 'destination' channel state
            is assigned with dst. A function that returns the transition rate in /s at 
            any voltage (in volts) is supplied with rate.
            
            Arguments: 
            * string id
            * steps.model.Surfsys surfsys
            * steps.model.ChanState src
            * steps.model.ChanState dst
            * function rate
            """
        ### Create 'ratelist' from input function
        if kwargs.has_key('vrange'):
            minv = kwargs['vrange'][0]
            maxv = kwargs['vrange'][1]
            dv = kwargs['vrange'][2]
            kwargs.pop('vrange')
        else:
            # These parameters are the hard-coded voltage bounds
            minv = -150.0e-3
            maxv = 100.0e-3
            dv = 1.0e-4
        
        tablesize = int((maxv-minv)/dv) +1
        ratelist = [0.0]*tablesize  
        rate = kwargs['rate']
        
        v = minv
        for i in range(tablesize):
        	ratelist[i] = rate(v)
        	v+=dv
        kwargs['ratetab'] = ratelist
        kwargs['vmin'] = minv
        kwargs['vmax'] = maxv
        kwargs['dv'] = dv
        kwargs['tablesize'] = tablesize
        kwargs.pop('rate')
        this = _steps_swig.new_VDepTrans(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        #self.__kwargs = kwargs
        self.__swig_setmethods__["id"] = _steps_swig.VDepTrans_setID
        self.__swig_getmethods__["id"] = _steps_swig.VDepTrans_getID
        self.__swig_getmethods__["model"] = _steps_swig.VDepTrans_getModel
        self.__swig_getmethods__["surfsys"] = _steps_swig.VDepTrans_getSurfsys        
        self.__swig_setmethods__["src"] = _steps_swig.VDepTrans_setSrc
        self.__swig_getmethods__["src"] = _steps_swig.VDepTrans_getSrc
        self.__swig_setmethods__["dst"] = _steps_swig.VDepTrans_setDst
        self.__swig_getmethods__["dst"] = _steps_swig.VDepTrans_getDst
    
        del(minv)
        del(maxv)
        del(dv)
        del(tablesize)
        del(ratelist)
        del(rate)
        del(v)
        
    id = steps_swig._swig_property(_steps_swig.VDepTrans_getID, _steps_swig.VDepTrans_setID)
    """Identifier string of the voltage-dependent transition."""
    model = steps_swig._swig_property(_steps_swig.VDepTrans_getModel)
    """Reference to parent model."""
    surfsys = steps_swig._swig_property(_steps_swig.VDepTrans_getSurfsys)
    """Reference to parent surface system."""
    src = steps_swig._swig_property(_steps_swig.VDepTrans_getSrc, _steps_swig.VDepTrans_setSrc)
    """ Reference to the channel state object that describes the 'source' channel. """ 
    dst = steps_swig._swig_property(_steps_swig.VDepTrans_getDst, _steps_swig.VDepTrans_setDst)
    """ Reference to the channel state object that describes the 'destination' channel. """ 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class VDepSReac(steps_swig.VDepSReac) :
    
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
    
    def __init__(self, *args, **kwargs): 
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
        ### Create 'ratelist' from input function
        if kwargs.has_key('vrange'):
            minv = kwargs['vrange'][0]
            maxv = kwargs['vrange'][1]
            dv = kwargs['vrange'][2]
            kwargs.pop('vrange')
        else:
            # These parameters are the hard-coded voltage bounds
            minv = -150.0e-3
            maxv = 100.0e-3
            dv = 1.0e-4
        
        tablesize = int((maxv-minv)/dv) +1
        klist = [0.0]*tablesize  
        k = kwargs['k']
        
        v = minv
        for i in range(tablesize):
        	klist[i] = k(v)
        	v+=dv
        kwargs['ktab'] = klist
        kwargs['vmin'] = minv
        kwargs['vmax'] = maxv
        kwargs['dv'] = dv
        kwargs['tablesize'] = tablesize
        kwargs.pop('k')
        this = _steps_swig.new_VDepSReac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        #self.__kwargs = kwargs
        self.__swig_setmethods__["id"] = _steps_swig.VDepSReac_setID
        self.__swig_getmethods__["id"] = _steps_swig.VDepSReac_getID
        self.__swig_getmethods__["model"] = _steps_swig.VDepSReac_getModel
        self.__swig_getmethods__["surfsys"] = _steps_swig.VDepSReac_getSurfsys        
        self.__swig_getmethods__["order"] = _steps_swig.VDepSReac_getOrder
        self.__swig_setmethods__["olhs"] = _steps_swig.VDepSReac_setOLHS
        self.__swig_getmethods__["olhs"] = _steps_swig.VDepSReac_getOLHS
        self.__swig_setmethods__["ilhs"] = _steps_swig.VDepSReac_setILHS
        self.__swig_getmethods__["ilhs"] = _steps_swig.VDepSReac_getILHS
        self.__swig_setmethods__["slhs"] = _steps_swig.VDepSReac_setSLHS
        self.__swig_getmethods__["slhs"] = _steps_swig.VDepSReac_getSLHS
        self.__swig_setmethods__["irhs"] = _steps_swig.VDepSReac_setIRHS
        self.__swig_getmethods__["irhs"] = _steps_swig.VDepSReac_getIRHS
        self.__swig_setmethods__["srhs"] = _steps_swig.VDepSReac_setSRHS
        self.__swig_getmethods__["srhs"] = _steps_swig.VDepSReac_getSRHS
        self.__swig_setmethods__["orhs"] = _steps_swig.VDepSReac_setORHS
        self.__swig_getmethods__["orhs"] = _steps_swig.VDepSReac_getORHS
    
        del(minv)
        del(maxv)
        del(dv)
        del(tablesize) 
        del(klist)  
        del(k)
        del(v)
    
    id = steps_swig._swig_property(_steps_swig.VDepSReac_getID, _steps_swig.VDepSReac_setID)
    """Identifier string of the voltage-dependent reaction."""
    model = steps_swig._swig_property(_steps_swig.VDepSReac_getModel)
    """Reference to parent model."""
    surfsys = steps_swig._swig_property(_steps_swig.VDepSReac_getSurfsys)
    """Reference to parent surface system."""
    olhs = steps_swig._swig_property(_steps_swig.VDepSReac_getOLHS, _steps_swig.VDepSReac_setOLHS)
    """Left hand side reactants in outer compartment."""
    ilhs = steps_swig._swig_property(_steps_swig.VDepSReac_getILHS, _steps_swig.VDepSReac_setILHS)
    """Left hand side reactants in inner compartment."""
    slhs = steps_swig._swig_property(_steps_swig.VDepSReac_getSLHS, _steps_swig.VDepSReac_setSLHS)
    """Left hand side reactants on surface."""
    irhs = steps_swig._swig_property(_steps_swig.VDepSReac_getIRHS, _steps_swig.VDepSReac_setIRHS)
    """Right hand side reactants in inner compartment."""
    srhs = steps_swig._swig_property(_steps_swig.VDepSReac_getSRHS, _steps_swig.VDepSReac_setSRHS)
    """Right hand side reactants on surface."""
    orhs = steps_swig._swig_property(_steps_swig.VDepSReac_getORHS, _steps_swig.VDepSReac_setORHS)
    """Right hand side reactants in outer compartment."""
    order = steps_swig._swig_property(_steps_swig.VDepSReac_getOrder)
    """Order of the voltage-dependent reaction."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class OhmicCurr(steps_swig.OhmicCurr) :
    
    """
        An ohmic current object that describes an ohmic current through a channel in 
        a particular conducting state. 
        
        An ohmic current object is a simple current based on a single, fixed value
        for single-channel conductance and a constant reversal potential. An ohmic current
        does not result in movement of ions between compartments, but simply 
        contributes a continous current to the EField solver for every conducting channel
        on a membrane surface that contains the ohmic current.
        """
    def __init__(self, *args, **kwargs): 
        """
            Construction::
            
            ohmiccurr = steps.model.OhmicCurr(id, surfsys, chanstate, erev, g) 
            
            Construct an ohmic curernt object with identifier string id
            and assign surfsys as the parent surface system. Assign to channel state
            chanstate, set the reversal potential to erev (in volts) and the single-channel
            conductance to g (in Siemens).
            
            Arguments: 
            * string id
            * steps.model.Surfsys surfsys
            * steps.model.ChanState chanstate
            * double erev
            * double g
            
            """
        this = _steps_swig.new_OhmicCurr(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.OhmicCurr_setID
        self.__swig_getmethods__["id"] = _steps_swig.OhmicCurr_getID
        self.__swig_getmethods__["model"] = _steps_swig.OhmicCurr_getModel
        self.__swig_getmethods__["surfsys"] = _steps_swig.OhmicCurr_getSurfsys        
        self.__swig_setmethods__["chanstate"] = _steps_swig.OhmicCurr_setChanState
        self.__swig_getmethods__["chanstate"] = _steps_swig.OhmicCurr_getChanState        
        self.__swig_setmethods__["erev"] = _steps_swig.OhmicCurr_setERev
        self.__swig_getmethods__["erev"] = _steps_swig.OhmicCurr_getERev
        self.__swig_setmethods__["g"] = _steps_swig.OhmicCurr_setG
        self.__swig_getmethods__["g"] = _steps_swig.OhmicCurr_getG
    
    id = steps_swig._swig_property(_steps_swig.OhmicCurr_getID, _steps_swig.OhmicCurr_setID)
    """Identifier string of the ohmic current."""
    model = steps_swig._swig_property(_steps_swig.OhmicCurr_getModel)
    """Reference to parent model."""
    surfsys = steps_swig._swig_property(_steps_swig.OhmicCurr_getSurfsys)
    """Reference to parent surface system."""
    erev = steps_swig._swig_property(_steps_swig.OhmicCurr_getERev, _steps_swig.OhmicCurr_setERev)
    """ The reversal potential (in volts). """
    g = steps_swig._swig_property(_steps_swig.OhmicCurr_getG, _steps_swig.OhmicCurr_setG)
    """ The single-channel conductance (in Siemens). """

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class GHKcurr(steps_swig.GHKcurr) :
    
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
    def __init__(self, *args, **kwargs): 
        """
            Construction::
            
            ghkcurr = steps.model.GHKcurr(id, surfsys, chanstate, ion, computeflux=True) 
            
            Construct a ghk current object with identifier string id
            and assign surfsys as the parent surface system. Assign to channel state
            chanstate, set the species that describes the current with ion- this
            species object must have a valence specified. If computeflux flag is 
            set to True then the current will result in movement of ions between compartments, 
            if False the current will be calculated but will not correspond to a real ion flux.
            
            Arguments: 
            * string id
            * steps.model.Surfsys surfsys
            * steps.model.ChanState chanstate
            * steps.model.Spec ion
            * bool computeflux
            
            NOTE: function setP or setPInfo must be called on the object before creating simulation object.
            """
        this = _steps_swig.new_GHKcurr(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.GHKcurr_setID
        self.__swig_getmethods__["id"] = _steps_swig.GHKcurr_getID
        self.__swig_getmethods__["model"] = _steps_swig.GHKcurr_getModel
        self.__swig_getmethods__["surfsys"] = _steps_swig.GHKcurr_getSurfsys        
        self.__swig_setmethods__["chanstate"] = _steps_swig.GHKcurr_setChanState
        self.__swig_getmethods__["chanstate"] = _steps_swig.GHKcurr_getChanState   
        self.__swig_setmethods__["ion"] = _steps_swig.GHKcurr_setIon
        self.__swig_getmethods__["ion"] = _steps_swig.GHKcurr_getIon      
        self.__swig_setmethods__["pinfo"] = _steps_swig.GHKcurr_setPInfo
        self.__swig_setmethods__["p"] = _steps_swig.GHKcurr_setP
    
    id = steps_swig._swig_property(_steps_swig.GHKcurr_getID, _steps_swig.GHKcurr_setID)
    """Identifier string of the ghk current."""
    model = steps_swig._swig_property(_steps_swig.GHKcurr_getModel)
    """Reference to parent model."""
    surfsys = steps_swig._swig_property(_steps_swig.GHKcurr_getSurfsys)
    """Reference to parent surface system."""
    ion = steps_swig._swig_property(_steps_swig.GHKcurr_getIon, _steps_swig.GHKcurr_setIon)
    """ The current ion. """ 
    pinfo = steps_swig._swig_property(_steps_swig.GHKcurr_setPInfo)       
    """ The infomation allowing permeability to be found internally. """
    p = steps_swig._swig_property(_steps_swig.GHKcurr_setP)       
    """ The permeability. """

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


