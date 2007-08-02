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
        self.__diffs = { }
        try:
            model._handleVolsysAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Volsys not assigned to model.'
            
            
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        for i in self.__reacs: 
            self.delReac(self.__reacs[i])
        self._model = None


    def _handleSpecDelete(self, spec):
        """
        """
        # Delete reaction rules using species.
        reacts = [ ]
        for i in self.__reacs: 
            r = self.__reacs[i]
            if (spec in r.lhs) or (spec in r.rhs):
                reacts.append(r)
        for r in reacts: self.delReac(r)
        # Delete diffusion rules for species.
        diffs = [ ]
        for i in self.__diffs:
            d = self.__difs[i]
            if spec in d.getAllSpecs():
                diffs.append(d)
        for d in diffs: self.delDiff(d)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _checkReacID(self, id):
        """Check if a given id is valid and not yet used by another reaction.

        PARAMETERS:
            id
                The id that should be checked (must be a string).

        RETURNS:
            The id itself.

        RAISES:
            steps.error.ArgumentError
                If the ID is not valid or unique.
        """
        stools.checkID(id)
        if id in self.__reacs:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id


    def _checkDiffID(self, id):
        """Check if a given id is valid and not yet used by another diffusion
        rule.
        
        PARAMETERS:
            id
                The id that should be checked (must be a string).
            
        RETURNS:
            The id itself.
        
        RAISES:
            steps.error.ArgumentError
                If the ID is not valid or unique.
        """
        stools.checkID(id)
        if id in self.__diffs:
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


    def getAllSpecs(self):
        """Create a list of all species involved in this volume system. The
        list contains no duplicates but is not ordered in any way.
        """
        s = set()
        for id, r in self.__reacs.iteritems():
            r = self.getReac(r)
            s = s.union(r.getAllSpecs())
        return list(s)
            
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleReacIDChange(self, oldid, newid):
        """
        """
        if oldid == newid: return
        self._checkReacID(newid)
        r = self.__reacs.pop(oldid)
        self.__reacs[newid] = r


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
                steps.volsys.Reac object.

        RETURNS:
            None

        RAISES:
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
                Can be a string or a steps.volsys.Reac object. If it is a
                string, the method will attempt to find the matching reaction
                object. If it is a steps.volsys.Reac object, the method 
                will see if the object is part of this volsys.

        RETURNS:
            A steps.volsys.Reac object.

        RAISES:
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
            raise serr.ArgumentError, "Reac is no part of this volsys."  
        return reaction
    

    def getAllReacs(self):
        """Return a list of all reactions defined in this volume system.
        """
        return self.__reacs.values()


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def _handleDiffIDChange(self, oldid, newid):
        """
        """
        if oldid == newid: return
        self._checkDiffID(newid)
        r = self.__diffs.pop(oldid)
        self.__diffs[newid] = r


    def _handleDiffAdd(self, diff):
        """
        """
        assert diff._volsys == None, \
            '\'%s\' already assigned to a volsys.' % diff.id
        self._checkDiffID(diff.id)
        diff._volsys = self
        self.__diffs[diff.id] = diff


    def delDiff(self, diff):
        """Remove a diffusion rule from the volume system.

        Arguments:
            diff
                Can be a string referring to the diffusion rule's id, or 
                a steps.volsys.Diff object.

        RETURNS:
            None

        RAISES:
            steps.error.ArgumentError
                If the id or diffusion rule object cannot be resolved.
        """
        diff = self.getDiff(diff)
        self.__diffs.pop(diff.id)
        diff._handleSelfDelete()


    def getDiff(self, diff):
        """Resolve a diffusion rule in the context of this volsys.

        This method can be used to retrieve a diffusion object by its
        id string, or to check whether a given diffusion object belongs
        to this volsys.

        Arguments:
            diff
                Can be a string or a steps.volsys.Diff object. If it is a
                string, the method will attempt to find the matching 
                diffusion object. If it is a steps.volsys.Diff object, 
                the method will see if the object is part of this volsys.

        RETURNS:
            A steps.volsys.Diff object.

        RAISES:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the diffusion rule 
                object does not belong to this model.
        """
        if isinstance(diff, str):
            try:
                diff = self.__diffs[diff]
            except KeyError:
                raise serr.ArgumentError, \
                    'No diffusion rule with id \'%s\'' % diff
        if diff._volsys != self:
            raise serr.ArgumentError, "Diff is no part of this volsys."  
        return diff
    

    def getAllDiffs(self):
        """Return a list of all diffusion rule object defined in this volsys.
        """
        return self.__diffs.values()
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Reac(object):


    """A reaction rule describing synthesis/degradation of chemical species.
    
    Currently, STEPS only supports the simplest form of reaction equations,
    in which reaction constants are really scalar constants: this means they 
    do not change .
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
        return len(self.__lhs)
    
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


class Diff(object):


    """A diffusion rule describing the diffusive motion of ligand species
    in a volume system.
    
    Currently, STEPS only supports the simplest form of diffusion, which
    is described by a single scalar constant. This means it doesn't depend
    on direction (i.e. is isotropic), doesn't change with time (unless 
    explicitly changed during simulation) or is not anomalous, and has no 
    dependencies on other simulation variables.
    
    In principle, it's possible to declare multiple diffusion rules for 
    one species in a volsys. This can be useful to simulate changing 
    diffusion properties, e.g. by activating only one diffusion rule for a 
    species at a time, and deactivating all others. (Admittedly, this is
    a dumb example -- or even a dumb feature -- right now because you might 
    as well just change the diffusion constant during simulation :-)
    """


    def __init__(self, id, volsys, lig, **params):
        """Initialize the object.
        """
        self._volsys = None
        self.__id = stools.checkID(id)
        try:
            volsys._handleDiffAdd(self)
        except:
            self._volsys = None
            raise
        assert self._volsys != None, 'Diffusion rule not assigned to volsys.'

        self.__lig = None
        self.lig = lig
        
        self._dcst = 0.0
        if 'dcst' in params: self.dcst = params['dcst']
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        self.__lig = None
        self._volsys = None

        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """Return ID of this diffusion rule.
        """
        return self.__id

    def setID(self, id):
        """
        """
        assert self._volsys != None, 'Diffusion rule not assigned to volsys.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._volsys._handleDiffIDChange(self.__id, id)
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
    
    
    def getLig(self):
        return self.__lig
    
    def setLig(self, lig):
        """Set the species to which this diffusion rule applies.
        
        This ligand can be changed after it has been specified in the
        initializer method.
        
        PARAMETERS:
            lig
                Can be a species object or a string id.
        """
        assert self._volsys != None
        lig = self.model.getSpec(lig)
        self.__lig = lig
    
    lig = property(getLig, setLig)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def getDcst(self):
        return self._dcst
    
    def setDcst(self, dcst):
        """Set the diffusion constant for the diffusion rule.
        
        PARAMETERS:
            dcst
                The diffusion constant, must be positive.
        
        TODO:
            Come up with a final solution for units of physical quantities.
        """
        dcst = float(dcst)
        if dcst < 0.0: 
            raise serr.ArgumentError, 'Diffusion constant must be >= 0.0'
        self._dcst = dcst

    dcst = property(getDcst, setDcst)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getAllSpecs(self):
        """Create a list of all species involved in this diffusion rule.
        
        Currently, this can obviously return only one species.
        """
        return [ self.__ligand ]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# END
