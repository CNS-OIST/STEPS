# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""
This file is the user-interface file for all rng objects in steps.
All objects are directly derived from the corresponding swig objects.
RNG object is owned by Python

"""
####### currently only imports swig file !! ##############

from .rng_swig import *

"""
import rng_swig
import _rng_swig

class RNG(rng_swig.RNG) :
    def __init__(self, *args, **kwargs): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _rng_swig.delete_RNG
    __del__ = lambda self : None;
    self.thisown = True

def create_mt19937(*args):
  return rng_swig.create_mt19937(*args)

def create(*args):
  return rng_swig.create(*args)
"""
