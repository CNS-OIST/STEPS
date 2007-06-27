# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""
"""


import steps.error as serr
import steps.tools as stools


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Volsys(object):
    

    """Container for reaction/diffusion rules occuring in volume solution.
    """
    

    def __init__(self, id, model):
        self._model = None
        self.__id = stools.checkID(id)
        self.__reacs = { }
        try:
            model._handleVolsysAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Volsys not assigned to model.'


    def _deepcopy(self, model):
        """
        """
        assert model != None
        v = steps.model.Volsys(self.id, model)
        for id, reac in self.__reacs.iteritems():
            reac._deepcopy(v)
            
            
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        for i, r in self.__reacs: self.delReac(r)
        self._model = None


    def _handleSpecDelete(self, spec):
        """
        """
        reacts = [ ]
        for i, r in self.__reacs.iteritems(): 
            if (spec in r.lhs) or (spec in r.rhs):
                reacts.append(r)
        for r in reacts: self.delReac(r)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _checkReacID(self, id):
        """Check if a given id is valid and not yet used by another reaction.

        Parameters:
            id
                The id that should be checked (must be a string).

        Returns:
            The id itself.

        Raises:
            steps.error.ArgumentError
                If the ID is not valid or unique.
        """
        stools.checkID(id)
        if id in self.__reacs:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """
        """
        return self.__id

    def setID(self, id):
        """
        """
        assert self._model != None, 'Volsys not assigned to model.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._model._handleVolsysIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getModel(self):
        """
        """
        return self._model

    model = property(getModel)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleReacIDChange(self, oldid, newid):
        """
        """
        if oldid == newid: return
        self._checkReacID(newid)
        r = self.__reacs.pop(oldid)
        self.__reacs[newid] = v


    def _handleReacAdd(self, reaction):
        """
        """
        assert reaction._volsys == None, \
            '\'%s\' already assigned to a volsys.' % reaction.id
        self._checkReacID(reaction.id)
        reaction._volsys = self
        self.__reacs[reaction.id] = reaction


    def delReac(self, reaction):
        """Remove a reaction from the volume system.

        Arguments:
            reaction
                Can be a string referring to the reaction's id, or a 
                steps.model.Reac object.

        Returns:
            None

        Raises:
            steps.error.ArgumentError
                If the id or reaction object cannot be resolved.
        """
        reaction = self.getReac(reaction)
        self.__reacs.pop(reaction.id)
        reaction._handleSelfDelete()


    def getReac(self, reaction):
        """Resolve a reaction rule in the context of this volsys.

        This method can be used to retrieve a reaction object by its
        id string, or to check whether a given reaction object belongs
        to this volsys.

        Arguments:
            reaction
                Can be a string or a steps.model.Reac object. If it is a
                string, the method will attempt to find the matching reaction
                object. If it is a steps.model.Reac object, the method 
                will see if the object is part of this volsys.

        Returns:
            A steps.model.Reac object.

        Raises:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the reaction object does
                not belong to this model.
        """
        if isinstance(reaction, str):
            try:
                reaction = self.__reacs[reaction]
            except KeyError:
                raise serr.ArgumentError, \
                    'No reaction with id \'%s\'' % reaction
        if reaction._volsys != self:
            raise steps.ModelError, "Reac is no part of this volsys."  
        return reaction
    

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getAllReacs(self):
        """Return a list of all reactions defined in this volume system.
        """
        return self.__reacs.values()


    def getAllSpecs(self):
        """Create a list of all species involved in this volume system. The
        list contains no duplicates but is not ordered in any way.
        """
        s = set()
        for id, r in self.__reacs.iteritems():
            r = self.getReac(r)
            s = s.union(r.getAllSpecs())
        return list(s)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Reac(object):


    """
    """

    
    def __init__(self, id, volsys, **params):
        """Initialize the object.
        """
        self._volsys = None
        self.__id = stools.checkID(id)
        try:
            volsys._handleReacAdd(self)
        except:
            self._volsys = None
            raise
        assert self._volsys != None, 'Reaction not assigned to volsys.'

        self.__lhs = [ ]
        if 'lhs' in params: self.lhs = params['lhs']
        self.__rhs = [ ]
        if 'rhs' in params: self.rhs = params['rhs']
        self._kcst = 0.0
        if 'kcst' in params: self.kcst = params['kcst']

    
    def _deepcopy(self, volsys):
        """
        """
        assert volsys != None
        r = steps.model.Reac(self.id, volsys)
        left = map(self.model.Spec.getID, self.__lhs)
        r.lhs = left
        right = map(self.model.Spec.getID, self.__rhs)
        r.rhs = right
        r.kcst = self.kcst
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        self.__lhs = [ ]
        self.__rhs = [ ]
        self._volsys = None

        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """Return ID of this reaction channel.
        """
        return self.__id

    def setID(self, id):
        """
        """
        assert self._volsys != None, 'Reaction not assigned to volsys.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._volsys._handleReacIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getVolsys(self):
        """
        """
        return self._volsys

    volsys = property(getVolsys)


    def getModel(self):
        """
        """
        return self._volsys.model

    model = property(getModel)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getLHS(self):
        """
        """
        assert self._volsys != None
        return self.__lhs[:]

    def setLHS(self, lhs):
        """
        """
        assert self._volsys != None
        # Turn into a list.
        if not isinstance(lhs, list): lhs = [lhs]
        lhs = map(self.model.getSpec, lhs)
        self.__lhs = lhs

    lhs = property(getLHS, setLHS)


    def getRHS(self):
        """
        """
        assert self._volsys != None
        return self.__rhs[:]

    def setRHS(self, rhs):
        """
        """
        assert self._volsys != None
        # Turn into list.
        if not isinstance(rhs, list): rhs = [rhs]
        rhs = map(self.model.getSpec, rhs)
        self.__rhs = rhs

    rhs = property(getRHS, setRHS)


    def getOrder(self):
        """
        """
        return len(lhs)
    
    order = property(getOrder)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getKcst(self):
        """
        """
        return self._kcst

    def setKcst(self, kcst):
        """
        """
        kcst = float(kcst)
        if kcst < 0.0: 
            raise serr.ArgumentError, 'Rate constant must be >= 0.0'
        self._kcst = kcst

    kcst = property(getKcst, setKcst)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getAllSpecs(self):
        """Create a list of all species involved in this reaction, on both
        the left- and righthand side. The list contains no duplicate members
        but is not ordered in any way.
        """
        s = set(self.__lhs)
        return list(s.union(self.__rhs))
        
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# END
