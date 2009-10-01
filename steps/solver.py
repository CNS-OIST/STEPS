# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2009 Stefan Wils. All rights reserved.
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
# $Id: solver.py 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""

This file is the user-interface file for all solver objects.
All objects are directly derived from the corresponding swig objects.
All objects are owned by Python.

"""

import solver_swig
import _solver_swig

class Wmrk4(solver_swig.Wmrk4) :
    def __init__(self, *args): 
        """__init__(self, steps::model::Model m, steps::wm::Geom g, steps::rng::RNG r) -> Wmrk4"""
        this = _solver_swig.new_Wmrk4(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True

class Wmdirect(solver_swig.Wmdirect) :
    def __init__(self, *args): 
        """__init__(self, steps::model::Model m, steps::wm::Geom g, steps::rng::RNG r) -> Wmdirect"""
        this = _solver_swig.new_Wmdirect(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True

class Tetexact(solver_swig.Tetexact) :
    def __init__(self, *args): 
        """__init__(self, steps::model::Model m, steps::wm::Geom g, steps::rng::RNG r) -> Tetexact"""
        this = _solver_swig.new_Tetexact(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True