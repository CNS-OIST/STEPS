# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Volsys(object):
    
    """Container for reaction/diffusion rules occuring in volume solution.
    """
    
    def __init__(self, id, model):
        self._model = None
        self.__id = steps.checkID(id)
        self.__reactions = { }
        try:
            model._handleVolsysAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Volsys not assigned to model.'

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleSelfDelete(self):
        for i, r in self.__reactions: self.delReaction(r)
        self._model = None

    def _handleSpeciesDelete(self, species):
        reacts = [ ]
        for i, r in self.__reactions.iteritems(): 
            if (species in r.lhs) or (species in r.rhs):
                reacts.append(r)
        for r in reacts: self.delReaction(r)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _checkReactionID(self, id):
        """Check if a given id is valid and not yet used by another reaction.

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
        if id in self.__reactions:
            raise steps.ModelError, '\'%s\' is already in use.' % id

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getID(self):
        return self.__id

    def setID(self, id):
        assert self._model != None, 'Volsys not assigned to model.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._model._handleVolsysIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getModel(self):
        return self._model

    model = property(getModel)
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleReactionIDChange(self, oldid, newid):
        if oldid == newid: return
        self._checkReactionID(newid)
        r = self.__reactions.pop(oldid)
        self.__reactions[newid] = v

    def _handleReactionAdd(self, reaction):
        assert reaction._volsys == None, \
            '\'%s\' already assigned to a volsys.' % reaction.id
        self._checkReactionID(reaction.id)
        reaction._volsys = self
        self.__reactions[reaction.id] = reaction

    def delReaction(self, reaction):
        """Remove a reaction from the volume system.

        Arguments:
            reaction
                Can be a string referring to the reaction's id, or a 
                steps.model.Reaction object.

        Returns:
            None

        Raises:
            steps.error.ModelError
                If the id or reaction object cannot be resolved.
        """
        reaction = self.getReaction(reaction)
        self.__reactions.pop(reaction.id)
        reaction._handleSelfDelete()

    def getReaction(self, reaction):
        """Resolve a reaction rule in the context of this volsys.

        This method can be used to retrieve a reaction object by its
        id string, or to check whether a given reaction object belongs
        to this volsys.

        Arguments:
            reaction
                Can be a string or a steps.model.Reaction object. If it is a
                string, the method will attempt to find the matching reaction
                object. If it is a steps.model.Reaction object, the method 
                will see if the object is part of this volsys.

        Returns:
            A steps.model.Reaction object.

        Raises:
            steps.error.ModelError
                If the id cannot be resolved, or if the reaction object does
                not belong to this model.
        """
        if isinstance(reaction, str):
            try:
                reaction = self.__reactions[reaction]
            except KeyError:
                raise steps.ModelError, \
                    'No reaction with id \'%s\'' % reaction
        if reaction._volsys != self:
            raise steps.ModelError, "Reaction is no part of this volsys."  
        return reaction
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getAllReactions(self):
        """Return a list of all reactions defined in this volume system."""
        return self.__reactions.values()

    def getAllSpecies(self):
        """Create a list of all species involved in this volume system. The
        list contains no duplicates but is not ordered in any way."""
        s = set()
        for id, r in self.__reactions.iteritems():
            r = self.getReaction(r)
            s = s.union(r.getAllSpecies())
        return list(s)
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def copy(self, model):
        assert model != None
        v = steps.model.Volsys(self.id, model)
        for id, reaction in self.__reactions.iteritems():
            reaction.copy(v)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
