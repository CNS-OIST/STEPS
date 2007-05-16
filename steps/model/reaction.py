# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Reaction(object):
    
    def __init__(self, id, volsys, **params):
        self._volsys = None
        self.__id = steps.checkID(id)
        try:
            volsys._handleReactionAdd(self)
        except:
            self._volsys = None
            raise
        assert self._volsys != None, 'Reaction not assigned to volsys.'

        self.__lhs = [ ]
        if 'lhs' in params: self.lhs = params['lhs']
        self.__rhs = [ ]
        if 'rhs' in params: self.rhs = params['rhs']
        self._unidir = False
        if 'unidir' in params: self.unidir = params['unidir']
        self._kf = 0.0
        if 'kf' in params: self.kf = params['kf']
        self._kb = 0.0
        if 'kb' in params: self.kb = params['kb']

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleSelfDelete(self):
        self.__lhs = [ ]
        self.__rhs = [ ]
        self._volsys = None
        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getID(self):
        return self.__id

    def setID(self, id):
        assert self._volsys != None, 'Reaction not assigned to volsys.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._volsys._handleReactionIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getVolsys(self):
        return self._volsys

    volsys = property(getVolsys)

    def getModel(self):
        return self._volsys.model

    model = property(getModel)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getLHS(self):
        assert self._volsys != None
        return self.__lhs[:]

    def setLHS(self, lhs):
        assert self._volsys != None
        # Turn into a list.
        if not isinstance(lhs, list): lhs = [lhs]
        try: lhs = map(self.model.getSpecies, lhs)
        except: pass
        else: self.__lhs = lhs

    lhs = property(getLHS, setLHS)

    def getForwardOrder(self):
        return len(lhs)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getRHS(self):
        assert self._volsys != None
        return self.__rhs[:]

    def setRHS(self, rhs):
        assert self._volsys != None
        # Turn into list.
        if not isinstance(rhs, list): rhs = [rhs]
        try: rhs = map(self.model.getSpecies, rhs)
        except: pass
        else: self.__rhs = rhs

    rhs = property(getRHS, setRHS)

    def getBackwardOrder(self):
        return len(rhs)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getUnidir(self):
        return self._unidir

    def setUnidir(self, ud):
        self._unidir = bool(ud)

    unidir = property(getUnidir, setUnidir)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getKf(self):
        return self._kf

    def setKf(self, kf):
        kf = float(kf)
        if kf < 0.0: raise ValueError, 'Rate constant must be >= 0.0'
        self._kf = kf

    kf = property(getKf, setKf)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getKb(self):
        return self._kb

    def setKb(self, kb):
        kb = float(kb)
        if kb < 0.0: raise ValueError, 'Rate constant must be >= 0.0'
        self._kb = kb

    kb = property(getKb, setKb)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getAllSpecies(self):
        """Create a list of all species involved in this reaction, on both
        the left- and righthand side. The list contains no duplicate members
        but is not ordered in any way."""
        s = set(self.__lhs)
        return list(s.union(self.__rhs))

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def copy(self, volsys):
        assert volsys != None
        r = steps.model.Reaction(self.id, volsys)
        left = map(self.model.Species.getID, self.__lhs)
        r.lhs = left
        right = map(self.model.Species.getID, self.__rhs)
        r.rhs = right
        r.unidir = self.unidir
        r.kf = self.kf
        r.kb = self.kb

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
