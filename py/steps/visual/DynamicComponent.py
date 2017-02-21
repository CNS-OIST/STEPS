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

from pyqtgraph.Qt import QtCore, QtGui, QtOpenGL
import pyqtgraph.opengl as gl
import numpy as np
import random

import steps.geom

class VisualTetsSpec(gl.GLScatterPlotItem):
    """
    Visualization component for species in tetrahedrons.
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * tets                List of tetrahedron indices
        * spec_id             ID of the species
        * spec_color          Color of the species
        * spec_size           Size of the species
        * max_nspec           Maximum number of species can be visualized in the component
        * max_density         Maximum density of species can be visualized in the component
        * auto_adjust         Boolean flag for auto adjustment of visualized species counts
    """
    def __init__(self, id, display, mesh, sim, tets, spec_id, spec_color = None, spec_size = 0.2, max_nspec = 10000, max_density = 1000e18, auto_adjust = True):
        """
        Constructor.
        """
        
        self.id = id
        self.display = display
        
        self.bound_min = [float('Inf'), float('Inf'), float('Inf')]
        self.bound_max = [-float('Inf'), -float('Inf'), -float('Inf')]
        
        self.mesh = mesh
        self.sim = sim
        self.spec_id = spec_id
        
        if spec_color == None:
            self.spec_color = [random.random(), random.random(), random.random(), random.random()]
        else:
            self.spec_color = spec_color
        
        self.spec_size = spec_size
        self.max_nspec = max_nspec
        self.max_density = max_density
        self.auto_adjust = auto_adjust
        
        self.tets = np.array(tets, dtype = np.uint32)
        self.counts = np.zeros(self.tets.size)
        sim.getBatchTetCountsNP(self.tets, spec_id, self.counts)
        
        point_counts = self.counts.astype(np.uint32)
        total = np.sum(point_counts)
        
        if total > max_nspec:
            total = self.__reduce(point_counts)
        
        data = np.zeros(total * 3)
        self.mesh.genTetVisualPointsNP(self.tets, point_counts, data)
        data *= display.scale
        
        data.shape = -1, 3
        
        gl.GLScatterPlotItem.__init__(self, pos = data, color=self.spec_color, size=self.spec_size, pxMode=False)
        display.addItem(self)
    
    def __reduce(self, point_counts):
        """
            Reduce the number of points being generated.
        """
        self.mesh.reduceBatchTetPointCountsNP(self.tets, point_counts, self.max_density)
        total = np.sum(point_counts)
        temp_density = self.max_density
        while total > self.max_nspec and self.auto_adjust:
            temp_density *= (self.max_nspec / total)
            self.mesh.reduceBatchTetPointCountsNP(self.tets, point_counts, temp_density)
            total = np.sum(point_counts)
        return total
    
    def updateItem(self):
        """
            Update the component.
        """
        self.sim.getBatchTetCountsNP(self.tets, self.spec_id, self.counts)
        point_counts = self.counts.astype(np.uint32)

        total = np.sum(point_counts)
        if total > self.max_nspec:
            total = self.__reduce(point_counts)
        data = np.zeros(total * 3)
        self.mesh.genTetVisualPointsNP(self.tets, point_counts, data)
        data *= self.display.scale
        data.shape = -1, 3
        self.setData(pos = data, color=self.spec_color)

    def setMaxNSpec(self, new_value):
        """
            Set the maximum number of species can be visualized in the component
        """
        self.max_nspec = new_value

    def getMaxNSpec(self):
        """
        Get the maximum number of species can be visualized in the component
        """
        return self.max_nspec

    def setSpecColor(self, color):
        """
        Set the color of the species
        """
        self.spec_color = color
        self.setData(color=self.spec_color)

    def getSpecColor(self):
        """
        Get the color of the species
        """
        return self.spec_color

    def setAutoAdjust(self, condition):
        """
        Set if the component automatically adjust the number of visualizing points according the maximum count and maximum density
        """
        self.auto_adjust = condition

    def getAutoAdjust(self):
        """
        Get if the component automatically adjust the number of visualizing points according the maximum count and maximum density
        """
        return self.auto_adjust

    def setMaxDensity(self, density):
        """
        Set the maximum density of species can be visualized in the component
        """
        self.max_density = density

    def getMaxDensity(self):
        """
        Get the maximum density of species can be visualized in the component
        """
        return self.max_density


    def setSpecSize(self, size):
        """
        Set the size of the points
        """
        self.spec_size = size
        self.setData(size = self.spec_size)

    def getSpecSize(self):
        """
        Get the size of the points
        """
        return self.spec_size

class VisualCompSpec(VisualTetsSpec):
    """
    Visualization component for species in compartment.
    
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * comp_id             ID of the compartment
        * spec_id             ID of the species
        * spec_color          Color of the species
        * spec_size           Size of the species
        * max_nspec           Maximum number of species can be visualized in the component
        * max_density         Maximum density of species can be visualized in the component
        * auto_adjust         Boolean flag for auto adjustment of visualized species counts
    """
    def __init__(self, id, display, mesh, sim, comp_id, spec_id, spec_color = None, spec_size = 0.2, max_nspec = 100000, max_density = 1000e18, auto_adjust = True):
        """
        Constructor.

        """
        self.comp_id = comp_id
        tets = steps.geom.castToTmComp(mesh.getComp(comp_id)).getAllTetIndices()
            
        VisualTetsSpec.__init__(self, id, display, mesh, sim, tets, spec_id, spec_color, spec_size, max_nspec, max_density, auto_adjust)

class VisualROITetsSpec(VisualTetsSpec):
    """
    Visualization component for species in Region of Interest tetrahedrons.
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * roi_id              ID of the Region of Interest
        * spec_id             ID of the species
        * spec_color          Color of the species
        * spec_size           Size of the species
        * max_nspec           Maximum number of species can be visualized in the component
        * max_density         Maximum density of species can be visualized in the component
        * auto_adjust         Boolean flag for auto adjustment of visualized species counts
    
    """
    def __init__(self, id, display, mesh, sim, roi_id, spec_id, spec_color = None, spec_size = 0.2, max_nspec = 100000, max_density = 1000e18, auto_adjust = True):
        """
        Constructor.
        """
        self.roi_id = roi_id
        tets = mesh.getROIData(roi_id)
        VisualTetsSpec.__init__(self, id, display, mesh, sim, tets, spec_id, spec_color, spec_size, max_nspec, max_density, auto_adjust)

class VisualTrisSpec(gl.GLScatterPlotItem):
    """
    Visualization component for species in triangles.
    
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * tris                List of triangle indices
        * spec_id             ID of the species
        * spec_color          Color of the species
        * spec_size           Size of the species
        * max_nspec           Maximum number of species can be visualized in the component
        * max_density         Maximum density of species can be visualized in the component
        * auto_adjust         Boolean flag for auto adjustment of visualized species counts
    """
    def __init__(self, id, display, mesh, sim, tris, spec_id, spec_color = None, spec_size = 0.2, max_nspec = 100000, max_density = 1000e12, auto_adjust = True):
        """
        Constructor.
        """
        self.id = id
        self.display = display
        
        self.bound_min = [float('Inf'), float('Inf'), float('Inf')]
        self.bound_max = [-float('Inf'), -float('Inf'), -float('Inf')]
        
        self.mesh = mesh
        self.sim = sim
        self.spec_id = spec_id
        
        if spec_color == None:
            self.spec_color = [random.random(), random.random(), random.random(), random.random()]
        else:
            self.spec_color = spec_color
        
        self.spec_size = spec_size
        self.max_nspec = max_nspec
        self.max_density = max_density
        self.auto_adjust = auto_adjust
        
        self.tris = np.array(tris, dtype = np.uint32)
        self.counts = np.zeros(self.tris.size)
        
        sim.getBatchTriCountsNP(self.tris, spec_id, self.counts)
        
        point_counts = self.counts.astype(np.uint32)
        total = np.sum(point_counts)
        if total > max_nspec:
            total = self.__reduce(point_counts)
        
        data = np.zeros(total * 3)
        self.mesh.genTriVisualPointsNP(self.tris, point_counts, data)
        data *= display.scale
        
        data.shape = -1, 3
        
        gl.GLScatterPlotItem.__init__(self, pos = data, color=self.spec_color, size=self.spec_size, pxMode=False)
        display.addItem(self)
    
    def __reduce(self, point_counts):
        """
            Reduce the number of points being generated.
        """
        self.mesh.reduceBatchTriPointCountsNP(self.tris, point_counts, self.max_density)
        total = np.sum(point_counts)
        temp_density = self.max_density
        while total > self.max_nspec and self.auto_adjust:
            temp_density *= (self.max_nspec / total)
            self.mesh.reduceBatchTriPointCountsNP(self.tris, point_counts, temp_density)
            total = np.sum(point_counts)
        return total
    
    def updateItem(self):
        """
            Update the component.
        """
        self.sim.getBatchTriCountsNP(self.tris, self.spec_id, self.counts)
        self.point_counts = self.counts.astype(np.uint32)
        self.total = np.sum(self.point_counts)
        if self.total > self.max_nspec:
            self.__reduce()
        
        data = np.zeros(self.total * 3)
        self.mesh.genTriVisualPointsNP(self.tris, self.point_counts, data)
        data *= self.display.scale
        data.shape = -1, 3
        self.setData(pos = data, color=self.spec_color)

    def setMaxNSpec(self, new_value):
        """
            Set the maximum number of species can be visualized in the component
        """
        self.max_nspec = new_value
            
    def getMaxNSpec(self):
        """
        Get the maximum number of species can be visualized in the component
        """
        return self.max_nspec
    
    def setSpecColor(self, color):
        """
        Set the color of the species
        """
        self.spec_color = color
        self.setData(color=self.spec_color)
    
    def getSpecColor(self):
        """
        Get the color of the species
        """
        return self.spec_color
    
    def setAutoAdjust(self, condition):
        """
        Set if the component automatically adjust the number of visualizing points according the maximum count and maximum density
        """
        self.auto_adjust = condition
    
    def getAutoAdjust(self):
        """
        Get if the component automatically adjust the number of visualizing points according the maximum count and maximum density
        """
        return self.auto_adjust
    
    def setMaxDensity(self, density):
        """
        Set the maximum density of species can be visualized in the component
        """
        self.max_density = density
    
    def getMaxDensity(self):
        """
        Get the maximum density of species can be visualized in the component
        """
        return self.max_density

    def setSpecSize(self, size):
        """
        Set the size of the points
        """
        self.spec_size = size
        self.setData(size = self.spec_size)
            
    def getSpecSize(self):
        """
        Get the size of the points
        """
        return self.spec_size

class VisualPatchSpec(VisualTrisSpec):
    """
    Visualization component for species in patch.
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * patch_id            ID of the patch
        * spec_id             ID of the species
        * spec_color          Color of the species
        * spec_size           Size of the species
        * max_nspec           Maximum number of species can be visualized in the component
        * max_density         Maximum density of species can be visualized in the component
        * auto_adjust         Boolean flag for auto adjustment of visualized species counts
    """
    def __init__(self, id, display, mesh, sim, patch_id, spec_id, spec_color = None, spec_size = 0.2, max_nspec = 100000, max_density = 1000e12, auto_adjust = True):

        """
        Constructor.

        """
        self.patch_id = patch_id
        tris = steps.geom.castToTmPatch(mesh.getPatch(patch_id)).getAllTriIndices()
        
        VisualTrisSpec.__init__(self, id, display, mesh, sim, tris, spec_id, spec_color, spec_size, max_nspec, max_density, auto_adjust)

class VisualROITrisSpec(VisualTrisSpec):
    """
    Visualization component for species in Region of Interest triangles.

    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * roi_id              ID of the Region of Interest
        * spec_id             ID of the species
        * spec_color          Color of the species
        * spec_size           Size of the species
        * max_nspec           Maximum number of species can be visualized in the component
        * max_density         Maximum density of species can be visualized in the component
        * auto_adjust         Boolean flag for auto adjustment of visualized species counts
    """
    def __init__(self, id, display, mesh, sim, roi_id, spec_id, spec_color = None, spec_size = 0.2, max_nspec = 100000, max_density = 1000e12, auto_adjust = True):
        """
        Constructor.
        """
        self.roi_id = roi_id
        tris = mesh.getROIData(roi_id)
        
        VisualTrisSpec.__init__(self, id, display, mesh, sim, tris, spec_id, spec_color, spec_size, max_nspec, max_density, auto_adjust)

class VisualTrisChannel(gl.GLScatterPlotItem):
    """
    Visualization component for channel species in triangles.
    
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * tris                List of triangle indices
        * specs_colors        Species-Color mapping dictionary
        * spec_size           Size of the species
    """
    def __init__(self, id, display, mesh, sim, tris, specs_colors, spec_size = 0.1):
        """
        Constructor.
        """
        self.id = id
        self.display = display
        
        self.bound_min = [float('Inf'), float('Inf'), float('Inf')]
        self.bound_max = [-float('Inf'), -float('Inf'), -float('Inf')]
        
        self.mesh = mesh
        self.sim = sim
        self.specs_colors = specs_colors
        
        self.spec_size = spec_size
        
        self.tris = np.array(tris, dtype = np.uint32)
        total_counts = np.zeros(self.tris.size)
        individual_counts = np.zeros(self.tris.size)
        
        for s in specs_colors.keys():
            self.sim.getBatchTriCountsNP(self.tris, s, individual_counts)
            total_counts += individual_counts
        
        self.starts = np.cumsum(total_counts) - 1
        total = sum(total_counts)
        data = np.zeros(total * 3)
        point_counts = total_counts.astype(np.uint32)
        self.mesh.genTriVisualPointsNP(self.tris, point_counts, data)
        data *= display.scale
        data.shape = -1, 3
        
        self.color = np.zeros(total * 4).reshape(-1, 4)
        
        current_counter = np.zeros(self.tris.size, dtype = np.uint32)
        for s in self.specs_colors.keys():
            sim.getBatchTriCountsNP(self.tris, s, individual_counts)
            for t in range(self.tris.size):
                self.color[self.starts[t] + current_counter[t] : self.starts[t] + current_counter[t] + individual_counts[t]] = self.specs_colors[s]
                current_counter[t] += individual_counts[t]

        
        gl.GLScatterPlotItem.__init__(self, pos = data, color=self.color, size=self.spec_size, pxMode=False)
        display.addItem(self)
    
    def updateItem(self):
        """
            Update the component.
        """
        current_counter = np.zeros(self.tris.size, dtype = np.uint32)
        individual_counts = np.zeros(self.tris.size)
        for s in self.specs_colors.keys():
            self.sim.getBatchTriCountsNP(self.tris, s, individual_counts)
            for t in range(self.tris.size):
                self.color[self.starts[t] + current_counter[t] : self.starts[t] + current_counter[t] + individual_counts[t]] = self.specs_colors[s]
                current_counter[t] += individual_counts[t]
        self.setData(color = self.color)
    
    def setSpecColor(self, spec, color):
        """
        Set the color of the species
        """
        self.specs_colors[spec] = color
        updateItem
    
    def getSpecColor(self, spec):
        """
        Get the color of the species
        """
        return self.specs_colors[spec]
    
    def setSpecSize(self, size):
        """
        Set the size of the points
        """
        self.spec_size = size
        self.setData(size = self.spec_size)
    
    def getSpecSize(self):
        """
        Get the size of the points
        """
        return self.spec_size

class VisualPatchChannel(VisualTrisChannel):
    """
    Visualization component for channel species in a patch.
    
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * patch_id            ID of the patch
        * specs_colors        Species-Color mapping dictionary
        * spec_size           Size of the species
    """
    def __init__(self, id, display, mesh, sim, patch_id, specs_colors, spec_size = 0.1):
        """
        Constructor.

        """
        self.patch_id = patch_id
        tris = steps.geom.castToTmPatch(mesh.getPatch(patch_id)).getAllTriIndices()
        
        VisualTrisChannel.__init__(self, id, display, mesh, sim, tris, specs_colors, spec_size)

class VisualROITrisChannel(VisualTrisChannel):
    """
    Visualization component for channel species in Region of Interest trangles.
    Parameters:
        * id                  ID of the component
        * display             Parent SimDisplay object
        * mesh                STEPS Tetmesh object
        * sim                 STEPS solver object
        * roi_id              ID of the Region of Interest
        * specs_colors        Species-Color mapping dictionary
        * spec_size           Size of the species
    """
    def __init__(self, id, display, mesh, sim, roi_id, specs_colors, spec_size = 0.1):
        """
        Constructor.
        """
        self.roi_id = roi_id
        tris = mesh.getROIData(roi_id)
        
        VisualTrisChannel.__init__(self, id, display, mesh, sim, tris, specs_colors, spec_size)

