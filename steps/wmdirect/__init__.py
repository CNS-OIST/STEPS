# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from steps.sim.controller import Controller
from wmdirect_core import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class WMDirect(Controller):

    """Controller class for Gillespie's Direct Method of SSA in well-mixed
    conditions.
    """

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def __init__(self, model, geom, rng):
        Controller.__init__(self, wmdirect_core, model, geom, rng)


    def __del__(self):
        Controller.__del__(self)
        

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
