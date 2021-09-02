####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pyqtgraph.opengl as gl
import time

QtCore.QCoreApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts)

class SimDisplay(QtGui.QMainWindow):
    """
    Simulation Display
    Parameters:

    * id                  ID of the display
    * x                   X cordinate of the display
    * y                   Y cordinate of the display
    * w                   Width of the display
    * h                   Height of the display
    * scale               Scaling between STEPS mesurement and display pixel
    """
    num_display = 0
    
    def __init__(self, id, x = 100, y = 100, w = 800, h = 600, scale = 1e6):
        """
        Constructor.
        """
        QtGui.QMainWindow.__init__(self)
        self.id = id
        self.setGeometry(x, y, w, h)
        self.widget = gl.GLViewWidget(self)
        self.setCentralWidget(self.widget)

        self.scale = scale
        self.center = [0.0,0.0,0.0]
        self.rotate = 0
        self.main_axis = 0
        
        self.setWindowTitle(id)
        
        self.bound_min = [float('Inf'), float('Inf'), float('Inf')]
        self.bound_max = [-float('Inf'), -float('Inf'), -float('Inf')]
        
        self.show()
    
    def closeEvent(self, event):
        self.hide()
        event.ignore()

    def addItem(self, item):
        """
        Add a visual component to the display.
            
        Parameters:

        * item                The adding visual component
            
        Return:
            None
        """
        self.widget.addItem(item)
        
        for i in range(3):
            if item.bound_min[i] < self.bound_min[i]:
                self.bound_min[i] = item.bound_min[i]
            if item.bound_max[i] > self.bound_max[i]:
                self.bound_max[i] = item.bound_max[i]

        self._updateView()
    
    def addItems(self, items):
        """
        Add a list of visual components to the display.
            
        Parameters:

        * items               List of the adding visual components
            
        Return:
            None
        """
        for item in items:
            self.widget.addItem(item)
        
            for i in range(3):
                if item.bound_min[i] < self.bound_min[i]:
                    self.bound_min[i] = item.bound_min[i]
                if item.bound_max[i] > self.bound_max[i]:
                    self.bound_max[i] = item.bound_max[i]

        self._updateView()
    
    def getItems(self):
        """
        Get a list of visual components in the display.
            
        Parameters:
            None
            
        Return:
            List of the visual components in the display
        """
        return self.widget.items
    
    def removeItem(self, item_id):
        """
        Remove a visual components in the display.
            
        Parameters:

        * item_id             ID of the removing component
            
        Return:
            None
        """
        for item in self.widget.items:
            if item.id == item_id:
                self.widget.items.remove(item)
                break
        
        self.bound_min = [float('Inf'), float('Inf'), float('Inf')]
        self.bound_max = [-float('Inf'), -float('Inf'), -float('Inf')]
        
        for item in self.widget.items:
            for i in range(3):
                if item.bound_min[i] < self.bound_min[i]:
                    self.bound_min[i] = item.bound_min[i]
                if item.bound_max[i] > self.bound_max[i]:
                    self.bound_max[i] = item.bound_max[i]
        
        self._updateView()
    
    
    def hideItem(self, item_id):
        """
        Hide a visual components in the display.
            
        Parameters:

        * item_id             ID of the hidding component
            
        Return:
            None
        """
        for item in self.widget.items:
            if item.id == item_id:
                item.hide()
    
    def showItem(self, item_id):
        """
        Show a visual components in the display.
            
        Parameters:

        * item_id             ID of the showing component
            
        Return:
            None
        """
        for item in self.widget.items:
            if item.id == item_id:
                item.show()
    
    def _updateView(self):
        
        axis = [0, 1, 2]
        new_center = [(self.bound_min[i] + self.bound_max[i]) / 2.0 for i in axis]
        
        pan_dist = [new_center[i] - self.center[i] for i in axis]
        self.widget.pan(pan_dist[0], pan_dist[1], pan_dist[2])
        
        dist_x = self.bound_max[0] - self.bound_min[0]
        dist_y = self.bound_max[1] - self.bound_min[1]
        dist_z = self.bound_max[2] - self.bound_min[2]
                
        if dist_y >= dist_x and dist_y >= dist_z:
            self.main_axis = 1
            self.widget.setCameraPosition(distance=dist_y * 2)
        
        elif dist_z >= dist_x and dist_z >= dist_y:
            self.widget.setCameraPosition(distance=dist_z * 2)
            self.main_axis = 2

        else:
            self.widget.setCameraPosition(distance=dist_x * 2)
            self.main_axis = 0

        self.center = new_center

    def updateItems(self):
        """
        Update visual components in the display.
            
        Parameters:
            None
            
        Return:
            None
        """
        for item in self.widget.items:
            if item.display == self:
                item.updateItem()
            else:
                item.update()

    def refresh(self):
        """
        Refresh the display.
            
        Parameters:
            None
            
        Return:
            None
        """
        self.widget.update()

    def rotateItems(self, angle, x, y, z, local=False):
        """
        Rotate all items in the display.
            
        Parameters:

        * angle           Angle for the rotation
        * x               X for the rotation vector
        * y               Y for the rotation vector
        * z               Z for the rotation vector
        * local           True if the rotation is local, False if it is global
            
        Return:
            None
        """
        for item in self.widget.items:
            item.rotate(angle, x, y, z, local)

    def rotateItem(self, item_id, angle, x, y, z, local=False):
        """
        Rotate all items in the display.
            
        Parameters:

        * item_id         ID of the rotating component
        * angle           Angle for the rotation
        * x               X for the rotation vector
        * y               Y for the rotation vector
        * z               Z for the rotation vector
        * local           True if the rotation is local, False if it is global
            
        Return:
            None
        """
        for item in self.widget.items:
            if item.id == item_id:
                item.rotate(angle, x, y, z, local)
                break

