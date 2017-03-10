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

__name__      = 'steps'
__longname__  = 'STochastic Engine for Pathway Simulation'
__version__   = '3.0.2'
__internal_version__ = 'HBP_0.10.8'
__author__    = 'STEPS Development Team'
__url__       = 'steps.sourceforge.net'
__license__   = 'GPL2.0'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

try:
    from . import steps_swig_numpy as steps_swig
    import _steps_swig_numpy as _steps_swig
except:
    print "Unable to import STEPS with NumPy support, try to import without NumPy..."
    try:
        from . import steps_swig
        import _steps_swig
    except:
        print "Unable to import STEPS."

_steps_swig.init()
import atexit
atexit.register(_steps_swig.finish)
import sys
if sys.stdout.isatty():
    _suppress_greet = False
    _quiet = False
else:
    _suppress_greet = True
    _quiet = True

def _greet():
    global _suppress_greet
    if not _suppress_greet:
        print ""
        print __longname__
        print "Version: ", __version__
        print "License: ", __license__
        print "Website: ", __url__

    _suppress_greet = True
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
