# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""Core classes and routines for describing model geometry.

A complete structural description of a model can be made with the classes 
in this module. Such descriptions can be used with well-mixed simulation
algorithms.
"""

from steps.tools import checkID

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Geom(object):

    def __init__(self):
        self.__comps = { }

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _checkCompID(self, id):
        checkID(id)
        if id in self.__comps:
            raise steps.GeomError, '\'%s\' is already in use.' % id

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleCompIDChange(self, oldid, newid):
        if oldid == newid: return
        self._checkCompID(newid)
        c = self.__comps.pop(oldid)
        self.__comps[newid] = c

    def _handleCompAdd(self, comp):
        assert comp._geom == None, \
            '\'%s\' already assigned to a geometry object.' % comp.id
        self._checkCompID(comp.id)
        comp._geom = self
        self.__comps[comp.id] = comp

    def delComp(self, comp):
        comp = self.getComp(comp)
        self.__comps.pop(comp.id)
        comp._handleSelfDelete()

    def getComp(self, comp):
        if isinstance(comp, str):
            try:
                comp = self.__comps[comp]
            except KeyError:
                raise steps.GeomError, 'No compartment with id \'%s\'' & comp
        if comp._geom != self:
            raise steps.GeomError, 'Compartment is no part of this geometry.'
        return comp
        
    def getAllComps(self):
        return self.__comps.values()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Comp(object):

    def __init__(self, id, geom, **params):
        self._geom = None
        self.__id = checkID(id)
        try:
            geom._handleCompAdd(self)
        except:
            self._geom = None
            raise
        assert self._geom != None, 'Compartment not assigned to geometry.'

        self.volsys = set()
        if 'volsys' in params: self.volsys = params['volsys']
        self._volume = 0.0
        if 'volume' in params: self._volume = params['volume']

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleSelfDelete(self):
        self._geom = None

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getID(self):
        return self.__id

    def setID(self, id):
        assert self._geom != None, 'Compartment not assigned to geometry.'
        if id == self.__id: return
        self._geom._handleCompIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getGeom(self):
        return self._geom

    geom = property(getGeom)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getVolume(self):
        return self._volume

    def setVolume(self, volume):
        # Python's float() type is actually a double..
        volume = float(volume)
        if volume <= 0.0:
            raise ValueError, 'Compartment volume must be larger than zero.'
        self._volume = volume

    volume = property(getVolume, setVolume)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
