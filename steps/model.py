# -*- coding: utf-8 -*-

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev$
#  Last Changed Date: $Date$
#  Last Changed By:   $Author$

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This file is the user-interface file for all model objects in steps.
# All objects are directly derived from the corresponding swig objects.
# Model container object is owned by Python
# All other objects are owned by c++ and Model container is responsible for 
# all the cleaning-up of these objects (see cpp/steps/model/model.cpp class 
# destructor).
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from . import model_swig
import _model_swig


###########################################
"""
class kcst(object):
    def __init__(self, fget, fset, _fget, _fset) :
        self.fget = fget
        self.fset = fset
        self._fget = _fget
        self._fset = _fset
    def __get__(self, instance, cls):
        k = self._fget(instance)
        self.fset(instance,k)
        return self.fget(instance)
    def __set__(self, instance, value):
        __swig_setmethods__(instance, value)
        self.fset(instance, value)
    def __swig_getmethods__ (self,instance,cls):
        return self._fget(instance)
    def __swig_setmethods__(self,instance,value):
        self._fset(instance,value)
"""
###########################################

class Model(model_swig.Model) : 
    """
    Top-level container for the objects in a kinetic model.
    
    """
    
    def __init__(self, *args):
        """
        Construction:
            m = steps.model.Model()
            
        Create a model container object.
            
        Arguments: 
            None
        """
        
        this = _model_swig.new_Model(*args)
        try: self.this.append(this)
        except: self.this = this
        # let the model object do the cleaning-up
        self.thisown = True

class Spec(model_swig.Spec) :

    """
    A chemical species which can be a reactant and/or product in reaction 
    stoichiometry and/or associated with a diffusion rule. 
    """
    def __init__(self, *args): 
        """        
        Construction:
            s = steps.model.Spec(id, mdl)
            
        Create a species object with identifier string id and assign the 
        object mdl as its parent model.
            
        Arguments: 
            * string id
            * steps.model.Model mdl
        """
        this = _model_swig.new_Spec(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Spec_setID
        self.__swig_getmethods__["id"] = _model_swig.Spec_getID
        self.__swig_getmethods__["model"] = _model_swig.Spec_getModel
    id = model_swig._swig_property(_model_swig.Spec_getID, _model_swig.Spec_setID)
    """Identifier string of the species."""
    model = model_swig._swig_property(_model_swig.Spec_getModel)
    """Reference to parent model."""

class Surfsys(model_swig.Surfsys):
    """
    A container that groups reactions involving a reactant embedded in a membrane.
    """
    def __init__(self, *args): 
        """
        Construction:
            s = steps.model.Surfsys(id, mdl)
            
        Construct a surface system object with identifier string id and assign 
        the object mdl as its parent model.
        
        Arguments: 
            * string id
            * steps.model.Model mdl
        """
        this = _model_swig.new_Surfsys(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Surfsys_setID
        self.__swig_getmethods__["id"] = _model_swig.Surfsys_getID
        self.__swig_getmethods__["model"] = _model_swig.Surfsys_getModel
    id = model_swig._swig_property(_model_swig.Surfsys_getID, _model_swig.Surfsys_setID)
    """Identifier string of the surface system."""
    model = model_swig._swig_property(_model_swig.Surfsys_getModel)
    """Reference to parent model."""

class Volsys(model_swig.Volsys) :
    """
    A container that groups reactions involving reactants embedded in a volume.
    """
    def __init__(self, *args): 
        """
        Construction:
            v = steps.model.Volsys(id, mdl)
            
        Construct a volume system object with identifier string id and assign 
        the object mdl as its parent model.
        
        Arguments: 
            * string id
            * steps.model.Model mdl
        """
        this = _model_swig.new_Volsys(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Volsys_setID
        self.__swig_getmethods__["id"] = _model_swig.Volsys_getID
        self.__swig_getmethods__["model"] = _model_swig.Volsys_getModel
    id = model_swig._swig_property(_model_swig.Volsys_getID, _model_swig.Volsys_setID)
    """Identifier string of the volume system."""
    model = model_swig._swig_property(_model_swig.Volsys_getModel)
    """Reference to parent model."""

class Diff(model_swig.Diff) :
    """
    A diffusion rule for a chemical species in a volume.
    
    A diffusion rule is described by:
    * Species to which the diffusion rule applies (lig).
    * Diffusion constant (dcst) specified in s.i. units.
    """
    def __init__(self, *args, **kwargs): 
        """
        Construction:
            diff = steps.model.Diff(id, volsys, lig, dcst = 0.0)
            
        Construct a diffusion rule object with identifier string id applied to species 
        lig and assign volsys as the parent volume system. Diffusion constant is 
        set by dcst.
        
        Arguments: 
            * string id
            * steps.model.Volsys volsys
            * steps.model.Spec lig
            * float dcst (default = 0.0)
        """
        this = _model_swig.new_Diff(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Diff_setID
        self.__swig_getmethods__["id"] = _model_swig.Diff_getID
        self.__swig_getmethods__["model"] = _model_swig.Diff_getModel
        self.__swig_getmethods__["volsys"] = _model_swig.Diff_getVolsys
        self.__swig_setmethods__["lig"] = _model_swig.Diff_setLig
        self.__swig_getmethods__["lig"] = _model_swig.Diff_getLig
        self.__swig_setmethods__["dcst"] = _model_swig.Diff_setDcst
        self.__swig_getmethods__["dcst"] = _model_swig.Diff_getDcst
    id = model_swig._swig_property(_model_swig.Diff_getID, _model_swig.Diff_setID)
    """Identifier string of the diffusion."""
    model = model_swig._swig_property(_model_swig.Diff_getModel)
    """Reference to parent model."""
    volsys = model_swig._swig_property(_model_swig.Diff_getVolsys)
    """Reference to parent volume system."""
    dcst = model_swig._swig_property(_model_swig.Diff_getDcst, _model_swig.Diff_setDcst)
    """Diffusion constant."""
    lig = model_swig._swig_property(_model_swig.Diff_getLig, _model_swig.Diff_setLig)
    """Reference to diffusion species."""

class Reac(model_swig.Reac) :
    """
    A kinetic reaction rule in a volume.
    
    The reaction rule is specified by:
    * Species on the left hand side of the reaction: the reactants (lhs).
    * Species on the right hand side of the reaction: the products (rhs).
    * Rate constant for the reaction, supplied in s.i. units (kcst).
    """
    def __init__(self, *args, **kwargs): 
        """
        Construction:
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
        this = _model_swig.new_Reac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Reac_setID
        self.__swig_getmethods__["id"] = _model_swig.Reac_getID
        self.__swig_getmethods__["model"] = _model_swig.Reac_getModel
        self.__swig_getmethods__["volsys"] = _model_swig.Reac_getVolsys
        self.__swig_setmethods__["kcst"] = _model_swig.Reac_setKcst
        self.__swig_getmethods__["kcst"] = _model_swig.Reac_getKcst
        self.__swig_setmethods__["lhs"] = _model_swig.Reac_setLHS
        self.__swig_getmethods__["lhs"] = _model_swig.Reac_getLHS
        self.__swig_setmethods__["rhs"] = _model_swig.Reac_setRHS
        self.__swig_getmethods__["rhs"] = _model_swig.Reac_getRHS
        self.__swig_getmethods__["order"] = _model_swig.Reac_getOrder
    id = model_swig._swig_property(_model_swig.Reac_getID, _model_swig.Reac_setID)
    """Identifier string of the reaction."""
    model = model_swig._swig_property(_model_swig.Reac_getModel)
    """Reference to parent model."""
    volsys = model_swig._swig_property(_model_swig.Reac_getVolsys)
    """Reference to parent volume system."""
    kcst = model_swig._swig_property(_model_swig.Reac_getKcst, _model_swig.Reac_setKcst)
    """Reaction constant."""
    lhs = model_swig._swig_property(_model_swig.Reac_getLHS, _model_swig.Reac_setLHS)
    """Left hand side reactants."""
    rhs = model_swig._swig_property(_model_swig.Reac_getRHS, _model_swig.Reac_setRHS)
    """Right hand side reactants."""
    order = model_swig._swig_property(_model_swig.Reac_getOrder)
    """Order of the reaction."""
    
class SReac(model_swig.SReac) :
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
        Construction:
            sreac = steps.model.Reac(id, surfsys, 
                                    ilhs = [ ], olhs = [ ], slhs = [ ],
                                    irhs = [ ], orhs = [], srhs = [ ],
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
        this = _model_swig.new_SReac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.SReac_setID
        self.__swig_getmethods__["id"] = _model_swig.SReac_getID
        self.__swig_getmethods__["model"] = _model_swig.SReac_getModel
        self.__swig_getmethods__["surfsys"] = _model_swig.SReac_getSurfsys
        self.__swig_setmethods__["kcst"] = _model_swig.SReac_setKcst
        self.__swig_getmethods__["kcst"] = _model_swig.SReac_getKcst
        self.__swig_getmethods__["order"] = _model_swig.SReac_getOrder
        self.__swig_getmethods__["outer"] = _model_swig.SReac_getOuter
        self.__swig_setmethods__["olhs"] = _model_swig.SReac_setOLHS
        self.__swig_getmethods__["olhs"] = _model_swig.SReac_getOLHS
        self.__swig_setmethods__["ilhs"] = _model_swig.SReac_setILHS
        self.__swig_getmethods__["ilhs"] = _model_swig.SReac_getILHS
        self.__swig_setmethods__["slhs"] = _model_swig.SReac_setSLHS
        self.__swig_getmethods__["slhs"] = _model_swig.SReac_getSLHS
        self.__swig_setmethods__["irhs"] = _model_swig.SReac_setIRHS
        self.__swig_getmethods__["irhs"] = _model_swig.SReac_getIRHS
        self.__swig_setmethods__["srhs"] = _model_swig.SReac_setSRHS
        self.__swig_getmethods__["srhs"] = _model_swig.SReac_getSRHS
        self.__swig_setmethods__["orhs"] = _model_swig.SReac_setORHS
        self.__swig_getmethods__["orhs"] = _model_swig.SReac_getORHS
    id = model_swig._swig_property(_model_swig.SReac_getID, _model_swig.SReac_setID)
    """Identifier string of the surface reaction."""
    model = model_swig._swig_property(_model_swig.SReac_getModel)
    """Reference to parent model."""
    surfsys = model_swig._swig_property(_model_swig.SReac_getSurfsys)
    """Reference to parent surface system."""
    kcst = model_swig._swig_property(_model_swig.SReac_getKcst, _model_swig.SReac_setKcst)
    """Reaction constant."""
    outer = model_swig._swig_property(_model_swig.SReac_getOuter)
    """Obsolete"""
    olhs = model_swig._swig_property(_model_swig.SReac_getOLHS, _model_swig.SReac_setOLHS)
    """Left hand side reactants in outer compartment."""
    ilhs = model_swig._swig_property(_model_swig.SReac_getILHS, _model_swig.SReac_setILHS)
    """Left hand side reactants in inner compartment."""
    slhs = model_swig._swig_property(_model_swig.SReac_getSLHS, _model_swig.SReac_setSLHS)
    """Left hand side reactants on surface."""
    irhs = model_swig._swig_property(_model_swig.SReac_getIRHS, _model_swig.SReac_setIRHS)
    """Right hand side reactants in inner compartment."""
    srhs = model_swig._swig_property(_model_swig.SReac_getSRHS, _model_swig.SReac_setSRHS)
    """Right hand side reactants on surface."""
    orhs = model_swig._swig_property(_model_swig.SReac_getORHS, _model_swig.SReac_setORHS)
    """Right hand side reactants in outer compartment."""
    order = model_swig._swig_property(_model_swig.SReac_getOrder)
    """Order of the reaction."""


