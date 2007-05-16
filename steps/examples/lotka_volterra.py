# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from steps.model import Model
from steps.model import Reaction
from steps.model import Species
from steps.model import Volsys

def make_model():
    lv = Model()
    vsys = Volsys('main', lv)
    

# END
