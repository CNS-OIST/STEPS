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


"""This module contains fundamental components for defining a STEPS model.

Currently, del() should not be called directly on the objects in this
module, because the garbage collector will not be able to collect them 
(due to references from the parent object). The only exception is the 
root object (in other words, objects of class Model).
"""


import steps.error as serr
import steps.tools as stools


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Model(object):


    """Top-level container for the objects in a kinetic model.
    
    A model.Model object is parent to the following objects:
        * model.Spec
        * volsys.Volsys
        * surfsys.Surfsys  
    """


    def __init__(self):
        """Initialize the Model object.
        """
        self.__specs = { }
        self.__volsys = { }
        self.__surfsys = { }
        
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    
    def _checkSpecID(self, id):
        """Check if a given id is valid and not yet used by another species.
        
        PARAMETERS:
            id
                The id that should be a checked (must be a string).

        RETURNS:
            The id itself.

        RAISES:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__specs:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id

    
    def _checkVolsysID(self, id):
        """Check if a given id is valid and not yet used by another volsys.

        PARAMETERS:
            id
                The id that should be checked (must be a string).

        RETURNS:
            The id itself.

        RAISES:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__volsys:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id


    def _checkSurfsysID(self, id):
        """Check if a given id is valid and not yet used by another surfsys.

        PARAMETERS:
            id
                The id that should be checked (must be a string).

        RETURNS:
            The id itself.

        RAISES:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__surfsys:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSpecIDChange(self, oldid, newid):
        """Handle a change in species ID.
        Internal method: called by method model.Spec.setID()
        """
        if oldid == newid: return
        self._checkSpecID(newid)
        s = self.__specs.pop(oldid)
        self.__specs[newid] = s


    def _handleSpecAdd(self, species):
        """Attempt to add a Spec object to this model.
        Internal method.
        """
        assert species._model == None, \
            '\'%s\' already assigned to a model.' % species.id
        self._checkSpecID(species.id)
        species._model = self
        self.__specs[species.id] = species


    def delSpec(self, species):
        """Remove a species from the entire model.

        Note: this will also remove all model components that refer to the
        species, such as reaction or diffusion rules.

        Arguments:
            species
                Can be a string referring to the id, or a 
                steps.model.Spec object.

        RETURNS:
            None

        RAISES:
            steps.error.ArgumentError
                If the id or species object cannot be resolved.
        """
        species = self.getSpec(species)
        for i in self.__volsys:
            v = self.__volsys[i] 
            v._handleSpecDelete(species)
        for i in self.__surfsys:
            s = self.__surfsys[i]
            v._handleSpecDelete(species)
        self.__specs.pop(species.id)
        species._handleSelfDelete()


    def getSpec(self, species):
        """Resolve a species in the context of this model.
        
        This method can be used to retrieve a species object by its id, or to
        check whether a given species object belongs to this model.
        
        Arguments:
            species
                Can be a string or a steps.model.Spec object. If it is a 
                string, the method will attempt to find the matching species 
                object. If it is a steps.model.Spec object, the method will
                see if the object is part of this model.

        RETURNS:
            A steps.model.Spec object.

        RAISES:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the species object does
                not belong to this model.
        """
        if isinstance(species, str):
            try:
                species = self.__specs[species]
            except KeyError:
                raise serr.ArgumentError, \
                    'No species with id \'%s\'' % species
        if species._model != self:
            raise serr.ArgumentError, 'Spec is no part of this model.'
        return species


    def getAllSpecs(self):
        """Return a list of all species.
        """
        return self.__specs.values()


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def _handleVolsysIDChange(self, oldid, newid):
        """Handle a change in volume system ID.
        Internal method.
        """
        if oldid == newid: return
        self._checkVolsysID(newid)
        v = self.__volsys.pop(oldid)
        self.__volsys[newid] = v


    def _handleVolsysAdd(self, volsys):
        """Handle the attempt to add a new volume system to the model.
        Internal method.
        """
        assert volsys._model == None, \
            '\'%s\' already assigned to a model.' % volsys.id
        self._checkVolsysID(volsys.id)
        volsys._model = self
        self.__volsys[volsys.id] = volsys


    def delVolsys(self, volsys):
        """Remove a volume system from the entire model.

        Arguments:
            volsys
                Can be a string referring to the id, or a 
                steps.model.Volsys object.

        RETURNS:
            None

        RAISES:
            steps.error.ModelError
                If the id or volsys object cannot be resolved.
        """
        volsys = self.getVolsys(volsys)
        self.__volsys.pop(volsys.id)
        volsys._handleSelfDelete()

        
    def getVolsys(self, volsys):
        """Resolve a volume system in the context of this model.
        
        This method can be used to retrieve a volsys object by its id, or to
        check whether a given volsys object belongs to this model.
        
        Arguments:
            volsys
                Can be a string or a steps.model.Volsys object. If it is a 
                string, the method will attempt to find the matching volsys 
                object. If it is a steps.model.Volsys object, the method will 
                see if the object is part of this model.

        RETURNS:
            A steps.model.Volsys object.

        RAISES:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the volsys object does
                not belong to this model.
        """
        if isinstance(volsys, str):
            try:
                volsys = self.__volsys[volsys]
            except KeyError:
                raise serr.ArgumentError, 'No volsys with id \'%s\'' % volsys
        if volsys._model != self:
            raise serr.ArgumentError, 'Volsys is no part of this model.'
        return volsys

        
    def getAllVolsys(self):
        """Return a list of all volume systems.
        """
        return self.__volsys.values()
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def _handleSurfsysIDChange(self, oldid, newid):
        """Handle a change in surface system ID.
        Internal method.
        """
        if oldid == newid: return
        self._checkSurfsysID(newid)
        s = self.__surfsys.pop(oldid)
        self.__surfsys[newid] = s


    def _handleSurfsysAdd(self, surfsys):
        """Handle the attempt to add a new surface system to the model.
        Internal method.
        """
        assert surfsys._model == None, \
            '\'%s\' already assigned to a model.' % surfsys.id
        self._checkSurfsysID(surfsys.id)
        surfsys._model = self
        self.__surfsys[surfsys.id] = surfsys


    def delSurfsys(self, surfsys):
        """Remove a surface system from the entire model.

        Arguments:
            surfsys
                Can be a string referring to the id, or a 
                steps.model.Surfsys object.

        RETURNS:
            None

        RAISES:
            steps.error.ModelError
                If the id or volsys object cannot be resolved.
        """
        surfsys = self.getSurfsys(surfsys)
        self.__surfsys.pop(surfsys.id)
        surfsys._handleSelfDelete()

        
    def getSurfsys(self, surfsys):
        """Resolve a surface system in the context of this model.
        
        This method can be used to retrieve a surfsys object by its id, or to
        check whether a given surfsys object belongs to this model.
        
        Arguments:
            surfsys
                Can be a string or a steps.model.Surfsys object. If it is a 
                string, the method will attempt to find the matching surfsys 
                object. If it is a steps.model.Surfsys object, the method will 
                see if the object is part of this model.

        RETURNS:
            A steps.model.Surfsys object.

        RAISES:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the surfsys object does
                not belong to this model.
        """
        if isinstance(surfsys, str):
            try:
                surfsys = self.__surfsys[surfsys]
            except KeyError:
                raise serr.ArgumentError, 'No surfsys with id \'%s\'' % surfsys
        if surfsys._model != self:
            raise serr.ArgumentError, 'Surfsys is no part of this model.'
        return surfsys

        
    def getAllSurfsys(self):
        """Return a list of all surface systems.
        """
        return self.__surfsys.values()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Spec(object):


    """A reactant that can be used in volume and surface systems. 
    """

    
    def __init__(self, id, model):
        """Initializer.
        """
        self._model = None
        self.__id = stools.checkID(id)
        try:
            model._handleSpecAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Species not assigned to model.'
        
        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """Handle a delete of this species.
        Internal method.
        """
        self._model = None


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """Return the species ID.
        """
        return self.__id

    def setID(self, id):
        """Set or change the species ID.
        """
        assert self._model != None, 'Species not assigned to model.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._model._handleSpecIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getModel(self):
        """Return a reference to the parent model.
        """
        return self._model

    model = property(getModel)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
