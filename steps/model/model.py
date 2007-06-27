# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""This module contains fundamental components for defining a STEPS model.
"""


import steps.error as serr
import steps.tools as stools


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Model(object):


    """Top-level container and management class for all model components.
    """


    def __init__(self):
        """Initialize the Model object.
        """
        self.__specs = { }
        self.__volsys = { }


    def deepcopy(self, visit):
        """Create a deep copy of the Model object, for use with standard
        module 'copy'.
        
        Deep copy means that not just the references, but all components
        belonging to the Model object themselves are copied recursively.
        """
        m = steps.model.Model()
        for id, spec in self.__specs.iteritems():
            spec._deepcopy(m)
        for id, volsys in self.__volsys.iteritems():
            volsys._deepcopy(m)
        return m
        
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    
    def _checkSpecID(self, id):
        """Check if a given id is valid and not yet used by another species.
        
        Parameters:
            id
                The id that should be a checked (must be a string).

        Returns:
            The id itself.

        Raises:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__specs:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id

    
    def _checkVolsysID(self, id):
        """Check if a given id is valid and not yet used by another volsys.

        Parameters:
            id
                The id that should be checked (must be a string).

        Returns:
            The id itself.

        Raises:
            steps.error.ArgumentError
                If the ID is not valid, or not unique.
        """
        stools.checkID(id)
        if id in self.__volsys:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSpecIDChange(self, oldid, newid):
        """
        """
        if oldid == newid: return
        self._checkSpecID(newid)
        s = self.__specs.pop(oldid)
        self.__specs[newid] = s


    def _handleSpecAdd(self, species):
        """
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

        Returns:
            None

        Raises:
            steps.error.ArgumentError
                If the id or species object cannot be resolved.
        """
        species = self.getSpec(species)
        for i, v in self.__volsys.iteritems(): 
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

        Returns:
            A steps.model.Spec object.

        Raises:
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
        """
        """
        return self.__specs.values()


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def _handleVolsysIDChange(self, oldid, newid):
        """
        """
        if oldid == newid: return
        self._checkVolsysID(newid)
        v = self.__volsys.pop(oldid)
        self.__volsys[newid] = v


    def _handleVolsysAdd(self, volsys):
        """
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

        Returns:
            None

        Raises:
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

        Returns:
            A steps.model.Volsys object.

        Raises:
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
        """
        """
        return self.__volsys.values()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Spec(object):


    """
    """

    
    def __init__(self, id, model):
        """
        """
        self._model = None
        self.__id = stools.checkID(id)
        try:
            model._handleSpecAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Species not assigned to model.'


    def _deepcopy(self, model):
        """
        """
        assert model != None
        s = steps.model.Spec(self.id, model)
        
        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        self._model = None


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """
        """
        return self.__id

    def setID(self, id):
        """
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
        """
        """
        return self._model

    model = property(getModel)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
