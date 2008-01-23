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


"""Surface systems describe the chemical environment on and near 2D
membranes. These will usually be real physical surfaces such as plasma 
membranes or internal structures.

During simulations, surface systems are linked to so-called patches (see
steps.geom package) in the same way that volume systems are linked to 
compartments. 
"""


import steps.error as serr
import steps.tools as stools


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Surfsys(object):
    

    """Container for kinetic rules occuring on surface membranes.
    """
    

    def __init__(self, id, model):
        self._model = None
        self.__id = stools.checkID(id)
        self.__sreacs = { }
        try:
            model._handleSurfsysAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Surfsys not assigned to model.'
            
            
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        for i in self.__ireacs: 
            self.delSReac(self.__sreacs[i])
        self._model = None


    def _handleSpecDelete(self, spec):
        """
        """
        reacts = [ ]
        for i in self.__sreacs: 
            r = self.__sreacs[i]
            if spec in r.getAllSpecs():
                reacts.append(r)
        for r in reacts: self.delSReac(r)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _checkSReacID(self, id):
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
        if id in self.__sreacs:
            raise serr.ArgumentError, '\'%s\' is already in use.' % id
        
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """
        """
        return self.__id

    def setID(self, id):
        """
        """
        assert self._model != None, 'Surfsys not assigned to model.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._model._handleSurfsysIDChange(self.__id, id)
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
        """Create a list of all species involved in this surface system. 
        The list contains no duplicates but is not ordered in any way.
        """
        s = set()
        for id, r in self.__sreacs.iteritems():
            r = self.getSReac(r)
            s = s.union(r.getAllSpecs())
        return list(s)
            
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSReacIDChange(self, oldid, newid):
        """
        """
        if oldid == newid: return
        self._checkSReacID(newid)
        r = self.__sreacs.pop(oldid)
        self.__sreacs[newid] = r


    def _handleSReacAdd(self, sreaction):
        """
        """
        assert sreaction._surfsys == None, \
            '\'%s\' already assigned to a surfsys.' % sreaction.id
        self._checkSReacID(sreaction.id)
        sreaction._surfsys = self
        self.__sreacs[sreaction.id] = sreaction


    def delSReac(self, sreaction):
        """Remove a surface reaction from the surface system.

        Arguments:
            sreaction
                Can be a string referring to the sreaction's id, or a 
                steps.surfsys.SReac object.

        RETURNS:
            None

        RAISES:
            steps.error.ArgumentError
                If the id or reaction object cannot be resolved.
        """
        sreaction = self.getSReac(sreaction)
        self.__sreacs.pop(sreaction.id)
        sreaction._handleSelfDelete()


    def getSReac(self, sreaction):
        """Resolve a sreaction rule in the context of this surfsys.

        This method can be used to retrieve a sreaction object by its
        id string, or to check whether a given sreaction object belongs
        to this surfsys.

        Arguments:
            sreaction
                Can be a string or a steps.surfsys.SReac object. If it is a
                string, the method will attempt to find the matching sreaction
                object. If it is a steps.surfsys.SReac object, the method 
                will see if the object is part of this surfsys.

        RETURNS:
            A steps.surfsys.SReac object.

        RAISES:
            steps.error.ArgumentError
                If the id cannot be resolved, or if the sreaction object does
                not belong to this model.
        """
        if isinstance(sreaction, str):
            try:
                sreaction = self.__sreacs[sreaction]
            except KeyError:
                raise serr.ArgumentError, \
                    'No sreaction with id \'%s\'' % sreaction
        if sreaction._surfsys != self:
            raise serr.ArgumentError, "SReac is no part of this surfsys."  
        return sreaction
    

    def getAllSReacs(self):
        """Return a list of all sreactions defined in this surface system.
        """
        return self.__sreacs.values()
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class SReac(object):
    
    
    """A reaction that involves one or more molecules embedded in surface
    membranes.
    
    The reaction can also involve volume reactants on the inside of the 
    membrane, or on the outside (but not on both sides). The products of
    a surface reaction can be on the outside volume, the inside volume or
    on the membrane itself (or any combination of this).
    """
    
    
    def __init__(self, id, surfsys, **params):
        """Initialize the object.
        """
        self._surfsys = None
        self.__id = stools.checkID(id)
        try:
            surfsys._handleSReacAdd(self)
        except:
            self._surfsys = None
            raise
        assert self._surfsys != None, \
            'Surface reaction not assigned to surfsys.'
        
        self.__outerlhs = True
        self.__vlhs = [ ]
        self.__slhs = [ ]
        self.__irhs = [ ]
        self.__srhs = [ ]
        self.__orhs = [ ]
        if 'vlhs' in params: self.vlhs = params['vlhs']
        if 'slhs' in params: self.slhs = params['slhs']
        if 'irhs' in params: self.irhs = params['irhs']
        if 'srhs' in params: self.srhs = params['srhs']
        if 'orhs' in params: self.orhs = params['orhs']
        self._kcst = 0.0
        if 'kcst' in params: self.kcst = params['kcst']
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def _handleSelfDelete(self):
        """
        """
        self.__vlhs = [ ]
        self.__slhs = [ ]
        self.__irhs = [ ]
        self.__srhs = [ ]
        self.__orhs = [ ]
        self._surfsys = None

        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getID(self):
        """Return ID of this surface reaction.
        """
        return self.__id

    def setID(self, id):
        """
        """
        assert self._surfsys != None, \
            'Surface reaction not assigned to surfsys.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._surfsys._handleSReacIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


    def getSurfsys(self):
        """
        """
        return self._surfsys

    surfsys = property(getSurfsys)


    def getModel(self):
        """
        """
        return self._surfsys.model

    model = property(getModel)


    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
    
    def getInner(self):
        return not self.__outerlhs
    
    def setInner(self, inner):
        self.__outerlhs = not bool(inner)
    
    inner = property(getInner, setInner)
    
    
    def getOuter(self):
        return self.__outerlhs
    
    def setOuter(self, outer):
        self.__outerlhs = bool(outer)
    
    outer = property(getOuter, setOuter)
    
    
    def getVLHS(self):
        """
        """
        assert self._surfsys != None
        return self.__vlhs[:]
    
    def setVLHS(self, vlhs):
        """
        """
        assert self._surfsys != None
        # Turn into a list.
        if not isinstance(vlhs, list): vlhs = [vlhs]
        vlhs = map(self.model.getSpec, vlhs)
        self.__vlhs = vlhs
    
    vlhs = property(getVLHS, setVLHS)
    
    
    def getSLHS(self):
        """
        """
        assert self._surfsys != None
        return self.__slhs[:]
    
    def setSLHS(self, slhs):
        """
        """
        assert self._surfsys != None
        # Turn into a list.
        if not isinstance(slhs, list): slhs = [slhs]
        slhs = map(self.model.getSpec, slhs)
        self.__slhs = slhs
    
    slhs = property(getSLHS, setSLHS)
    
    
    def getIRHS(self):
        """
        """
        assert self._surfsys != None
        return self.__irhs[:]
    
    def setIRHS(self, irhs):
        """
        """
        assert self._surfsys != None
        # Turn into list.
        if not isinstance(irhs, list): irhs = [irhs]
        irhs = map(self.model.getSpec, irhs)
        self.__irhs = irhs
    
    irhs = property(getIRHS, setIRHS)
    
    
    def getSRHS(self):
        """
        """
        assert self._surfsys != None
        return self.__srhs[:]
    
    def setSRHS(self, srhs):
        """
        """
        assert self._surfsys != None
        # Turn into list.
        if not isinstance(srhs, list): srhs = [srhs]
        srhs = map(self.model.getSpec, srhs)
        self.__srhs = srhs
    
    srhs = property(getSRHS, setSRHS)
    
    
    def getORHS(self):
        """
        """
        assert self._surfsys != None
        return self.__orhs[:]
    
    def setORHS(self, orhs):
        """
        """
        assert self._surfsys != None
        # Turn into list.
        if not isinstance(orhs, list): orhs = [orhs]
        orhs = map(self.model.getSpec, orhs)
        self.__orhs = orhs
    
    orhs = property(getORHS, setORHS)
    
    
    def getOrder(self):
        """
        """
        return len(self.__vlhs) + len(self.__slhs)
    
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
        """Create a list of all species involved in this sreaction, on both
        the left- and righthand side. The list contains no duplicate members
        but is not ordered in any way.
        """
        s = set(self.__vlhs)
        s.union(self.__slhs)
        s.union(self.__irhs)
        s.union(self.__srhs)
        s.union(self.__orhs)
        return list(s)        


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# END
