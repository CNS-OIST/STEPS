# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from steps.sim.controller import FuncCore, FuncSSA, FuncTetmesh
import steps.sim.tetexact_core as tetexact_core

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class TetExact(FuncTetmesh, FuncSSA, FuncCore):

    """Controller class for Gillespie's Direct Method of SSA over a
    tetrahedral mesh.
    """

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def __init__(self, model, geom, rng):
        FuncCore.__init__(self, tetexact_core, model, geom, rng)
        FuncSSA.__init__(self)
        FuncTetmesh.__init__(self, geom)


    def __del__(self):
        FuncSSA.__del__(self)
        FuncCore.__del__(self)
        FuncTetmesh.__del__(self)
        

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
