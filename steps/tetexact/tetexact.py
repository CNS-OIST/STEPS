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
