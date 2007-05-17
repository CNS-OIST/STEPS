# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from steps.sim.controller import FuncCore, FuncSSA
import steps.sim.wmdirect_core as wmdirect_core

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class WMDirect(FuncSSA, FuncCore):

    """Controller class for Gillespie's Direct Method of SSA in well-mixed
    conditions.
    """

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def __init__(self, model, geom, rng):
        FuncCore.__init__(self, wmdirect_core, model, geom, rng)
        FuncSSA.__init__(self)


    def __del__(self):
        FuncSSA.__del__(self)
        FuncCore.__del__(self)
        

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
