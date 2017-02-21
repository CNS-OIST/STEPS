# -*- coding: utf-8 -*-

####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

try:
    from . import steps_swig_numpy as steps_swig
    import _steps_swig_numpy as _steps_swig
except:
    from . import steps_swig
    import _steps_swig

import steps

### Now defunct mesh saving/loading tool ###
# from steps_swig import loadASCII, saveASCII


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class RNG(steps_swig.RNG):
    """
        Base class for all random number generators in STEPS.
    """
    def __init__(self, *args): 
        """
        """
        this = _steps_swig.new_RNG(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def create(*args):
  """
    Creates and returns a reference to a steps.rng.RNG random number generator object, 
    which is specified by type and pre-allocates a buffer list with size of buffer_size.

    Syntax::
        
        create(type, buffer_size)

    Arguments:
        * string type
        * uint buffer_size

    Return:
        steps.rng.RNG

    """
  return steps_swig.create(*args)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def create_mt19937(*args):
  """
    Creates and returns a reference to a steps.rng.RNG random number generator object, 
    which is specified by type and pre-allocates a buffer list with size of buffer_size.

    Syntax::
        
        create_mt19937(buffer_size)

    Arguments:
        * uint buffer_size

    Return:
        steps.rng.RNG

    """
  return steps_swig.create_mt19937(*args)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
