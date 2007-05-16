# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Species(object):
    
    def __init__(self, id, model):
        self._model = None
        self.__id = steps.checkID(id)
        try:
            model._handleSpeciesAdd(self)
        except:
            self._model = None
            raise
        assert self._model != None, 'Species not assigned to model.'

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def _handleSelfDelete(self):
        self._model = None

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getID(self):
        return self.__id

    def setID(self, id):
        assert self._model != None, 'Species not assigned to model.'
        if id == self.__id: return
        # The following might raise an exception, e.g. if the id is not
        # valid or not unique.
        self._model._handleSpeciesIDChange(self.__id, id)
        self.__id = id

    id = property(getID, setID)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def getModel(self):
        return self._model

    model = property(getModel)

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

    def copy(self, model):
        assert model != None
        s = steps.model.Species(self.id, model)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
