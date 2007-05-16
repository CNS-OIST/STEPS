# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Model(object):

    """The top-level container and management class for all model components.
    """

    def __init__(self):
        self.__species = { }
        self.__volsys = { }
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    def _checkSpeciesID(self, id):
        """Check if a given id is valid and not yet used by another species.
        
        Parameters:
            id
                The id that should be a checked (must be a string).

        Returns:
            The id itself.

        Raises:
            steps.error.IDError
                If the ID is not valid.
            steps.error.ModelError
                If the ID is not unique.
        """
        steps.checkID(id)
        if id in self.__species:
            raise steps.ModelError, '\'%s\' is already in use.' % id
    
    def _checkVolsysID(self, id):
        """Check if a given id is valid and not yet used by another volsys.

        Parameters:
            id
                The id that should be checked (must be a string).

        Returns:
            The id itself.

        Raises:
            steps.error.IDError
                If the ID is not valid.
            steps.error.ModelError
                If the ID is not unique.
        """
        steps.checkID(id)
        if id in self.__volsys:
            raise steps.ModelError, '\'%s\' is already in use.' % id

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleSpeciesIDChange(self, oldid, newid):
        if oldid == newid: return
        self._checkSpeciesID(newid)
        s = self.__species.pop(oldid)
        self.__species[newid] = s

    def _handleSpeciesAdd(self, species):
        assert species._model == None, \
            '\'%s\' already assigned to a model.' % species.id
        self._checkSpeciesID(species.id)
        species._model = self
        self.__species[species.id] = species

    def delSpecies(self, species):
        """Remove a species from the entire model.

        Note: this will also remove all model components that refer to the
        species, such as reaction or diffusion rules.

        Arguments:
            species
                Can be a string referring to the id, or a 
                steps.model.Species object.

        Returns:
            None

        Raises:
            steps.error.ModelError
                If the id or species object cannot be resolved.
        """
        species = self.getSpecies(species)
        for i, v in self.__volsys.iteritems(): 
            v._handleSpeciesDelete(species)
        self.__species.pop(species.id)
        species._handleSelfDelete()

    def getSpecies(self, species):
        """Resolve a species in the context of this model.
        
        This method can be used to retrieve a species object by its id, or to
        check whether a given species object belongs to this model.
        
        Arguments:
            species
                Can be a string or a steps.model.Species object. If it is a 
                string, the method will attempt to find the matching species 
                object. If it is a steps.model.Species object, the method will
                see if the object is part of this model.

        Returns:
            A steps.model.Species object.

        Raises:
            steps.error.ModelError
                If the id cannot be resolved, or if the species object does
                not belong to this model.
        """
        if isinstance(species, str):
            try:
                species = self.__species[species]
            except KeyError:
                raise steps.ModelError, 'No species with id \'%s\'' % species
        if species._model != self:
            raise steps.ModelError, 'Species is no part of this model.'
        return species

    def getAllSpecies(self):
        return self.__species.values()

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    def _handleVolsysIDChange(self, oldid, newid):
        if oldid == newid: return
        self._checkVolsysID(newid)
        v = self.__volsys.pop(oldid)
        self.__volsys[newid] = v

    def _handleVolsysAdd(self, volsys):
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
            steps.error.ModelError
                If the id cannot be resolved, or if the volsys object does
                not belong to this model.
        """
        if isinstance(volsys, str):
            try:
                volsys = self.__volsys[volsys]
            except KeyError:
                raise steps.ModelError, 'No volsys with id \'%s\'' % volsys
        if volsys._model != self:
            raise steps.ModelError, 'Volsys is no part of this model.'
        return volsys
        
    def getAllVolsys(self):
        return self.__volsys.values()
        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def copy(self):
        m = steps.model.Model()
        for id, species in self.__species:
            species.copy(m)
        for id, volsys in self.__volsys:
            volsys.copy(m)
        return m

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
