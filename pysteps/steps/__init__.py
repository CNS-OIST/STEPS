# -*- coding: utf-8 -*-

####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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
from __future__ import print_function, absolute_import

#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import glob
import atexit
import os.path

__name__      = 'steps'
__longname__  = 'STochastic Engine for Pathway Simulation'
__version__   = '3.3.0'
__author__    = 'STEPS Development Team'
__url__       = 'steps.sourceforge.net'
__license__   = 'GPL2.0'
__binding__   = 'Cython'

#Try importing mpi version if file exists. Avoids catching exception for the wrong reason
if glob.glob(os.path.join(os.path.dirname(__file__),'cysteps_mpi.*')):
    from . import cysteps_mpi
    from . import cysteps_mpi as stepslib
elif glob.glob(os.path.join(os.path.dirname(__file__),'cysteps.*')):
    from . import cysteps
    from . import cysteps as stepslib
else:
    raise Exception("STEPS library not found [cysteps*]")

stepslib._py_init()
atexit.register(stepslib._py_finish)

_suppress_greet = False
_quiet = False

def _greet():
    global _suppress_greet
    if not _suppress_greet:
        print("")
        print(__longname__)
        print("Version: ", __version__)
        print("License: ", __license__)
        print("Website: ", __url__)
        print("CXX Binding:", __binding__)

    _suppress_greet = True
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
