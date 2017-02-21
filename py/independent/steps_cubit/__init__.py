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

__name__      = 'steps_cubit'
__longname__  = 'CUBIT Mesh Preparation Support Toolkit for STEPS'
__version__   = '1.0.0'
__author__    = 'STEPS Development Team'
__url__       = 'steps.sourceforge.net'
__license__   = 'GPL2.0'
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

print ""
print ""
print "####################", __longname__, "#################"
print "Version: ", __version__
print "License: ", __license__
print "Website: ",  __url__
print ""
print ""

try:
    from steps_cubit import *
except:
    pass
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
