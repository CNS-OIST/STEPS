####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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
from __future__ import print_function
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from pyqtgraph.Qt import QtCore, QtGui

def pickColorF():
    """
    Pick a color and return it as a (0,1) floating point [r,g,b,a] list
        
    Parameters:
        None
        
    Return:
        List of the color
    """
    col = QtGui.QColorDialog.getRgba()
    qcolor = QtGui.QColor.fromRgba(col[0])
    return [qcolor.redF(), qcolor.greenF(), qcolor.blueF(), qcolor.alphaF()]

def printColorF():
    """
    Pick a color and print(out the list)
        
    Parameters:
        None
        
    Return:
        None
    """
    print(pickColorF())
