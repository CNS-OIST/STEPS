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

import pyqtgraph as pg
import numpy as np
import steps.geom

class PlotDisplay(pg.GraphicsWindow):
    """
    Visualization plot display.
    
    Parameters:
        * title           Title of the display
        * size            Size of the display
    """
    def __init__(self, title=None, size=(800, 600)):
        """
        Constructor.
        """
        pg.GraphicsWindow.__init__(self, title=title, size=size)
        self.updater = {}
    
    def addCompSpecPlot(self, title, sim, comp_id, spec_id, data_size = 1000, x_range = None, y_range = None, measure = "count", **kwargs):
        """
        Add plot to display the amount changes of species in a compartment.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * comp_id         ID of the compartment
            * spec_id         ID of the species
            * data_size       Size of the data history
            * x_range         Range of X axis
            * y_range         Range of y axis
            * measure         Measure type
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = CompSpecUpdater(plot, sim, comp_id, spec_id, data_size, x_range, y_range, measure, **kwargs)
        self.updater[title] = updater
        return plot
        
    def addTetsSpecPlot(self, title, sim, tets, spec_id, data_size = 1000, x_range = None, y_range = None, measure = "count", **kwargs):
        """
        Add plot to display the amount changes of species in a list of tetrahedrons.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * tets            List of tetrahedron indices
            * spec_id         ID of the species
            * data_size       Size of the data history
            * x_range         Range of X axis
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = TetsSpecUpdater(plot, sim, tets, spec_id, data_size, x_range, y_range, measure, **kwargs)
        self.updater[title] = updater
        return plot
    
    def addPatchSpecPlot(self, title, sim, patch_id, spec_id, data_size = 1000, x_range = None, y_range = None, **kwargs):
        """
        Add plot to display the amount changes of species in a patch.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * patch_id        ID of the patch
            * spec_id         ID of the species
            * data_size       Size of the data history
            * x_range         Range of X axis
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = PatchSpecUpdater(plot, sim, patch_id, spec_id, data_size, x_range, y_range, **kwargs)
        self.updater[title] = updater
        return plot

    def addTrisSpecPlot(self, title, sim, tris, spec_id, data_size = 1000, x_range = None, y_range = None, **kwargs):
        """
        Add plot to display the amount changes of species in a list of triangles.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * tris            List of triangle indices
            * spec_id         ID of the species
            * data_size       Size of the data history
            * x_range         Range of X axis
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = TrisSpecUpdater(plot, sim, tris, spec_id, data_size, x_range, y_range, **kwargs)
        self.updater[title] = updater
        return plot

    def addCompSumSpecsPlot(self, title, sim, comp_id, spec_ids, data_size = 1000, x_range = None, y_range = None, **kwargs):
        """
        Add plot to display the sum-up amount changes of species in a compartment.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * comp_id         ID of the compartment
            * spec_id         ID of the species
            * data_size       Size of the data history
            * x_range         Range of X axis
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = CompSumSpecsUpdater(plot, sim, comp_id, spec_ids, data_size, x_range, y_range, **kwargs)
        self.updater[title] = updater
        return plot
    
    def addPatchSumSpecsPlot(self, title, sim, patch_id, spec_ids, data_size = 1000, x_range = None, y_range = None, **kwargs):
        """
        Add plot to display the sum-up amount changes of species in a patch.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * patch_id        ID of the patch
            * spec_id         ID of the species
            * data_size       Size of the data history
            * x_range         Range of X axis
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = PatchSumSpecsUpdater(plot, sim, patch_id, spec_ids, data_size, x_range, y_range, **kwargs)
        self.updater[title] = updater
        return plot
    
    def addCompSpecDist(self, title, mesh, sim, comp_id, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        """
        Add plot to display the distribution of species in a compartment.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * comp_id         ID of the compartment
            * spec_id         ID of the species
            * axis            Spatial direction of the distribution
            * nbins           Number of bins for the data
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = CompSpecDistUpdater(plot, mesh, sim, comp_id, spec_id, axis, nbins, y_range, **kwargs)
        self.updater[title] = updater
        return plot
        
    def addPatchSpecDist(self, title, mesh, sim, patch_id, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        """
        Add plot to display the distribution of species in a patch.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * mesh            STEPS mesh
            * patch_id        ID of the patch
            * spec_id         ID of the species
            * axis            Spatial direction of the distribution
            * nbins           Number of bins for the data
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = PatchSpecDistUpdater(plot, mesh, sim, patch_id, spec_id, axis, nbins, y_range, **kwargs)
        self.updater[title] = updater
        return plot
    
    def addTetsSpecDist(self, title, mesh, sim, tets, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        """
        Add plot to display the distribution of species in a list of tetrahedrons.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * mesh            STEPS mesh
            * tets            List of tetrahedron indices
            * spec_id         ID of the species
            * axis            Spatial direction of the distribution
            * nbins           Number of bins for the data
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = TetsSpecDistUpdater(plot, mesh, sim, tets, spec_id, axis, nbins, y_range, **kwargs)
        self.updater[title] = updater
        return plot
        
    def addTrisSpecDist(self, title, mesh, sim, tris, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        """
        Add plot to display the distribution of species in a list of triangles.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * mesh            STEPS mesh
            * tris            List of triangle indices
            * spec_id         ID of the species
            * axis            Spatial direction of the distribution
            * nbins           Number of bins for the data
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        updater = TrisSpecDistUpdater(plot, mesh, sim, tris, spec_id, axis, nbins, y_range, **kwargs)
        self.updater[title] = updater
        return plot
        
    def addROISpecDist(self, title, mesh, sim, roi_id, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        """
        Add plot to display the distribution of species in a Region of Interest.
            
        Parameters:
            * title           Title of the plot
            * sim             STEPS solver
            * mesh            STEPS mesh
            * roi_id          ID of the Region of Interest
            * spec_id         ID of the species
            * axis            Spatial direction of the distribution
            * nbins           Number of bins for the data
            * y_range         Range of y axis
            * **kwargs        Other keywords that are supported by pygraph.GraphicsWindow class
            
        Return:
            pyqtgraph.PlotItem object
        """
        if title in self.updater:
            raise NameError('A Plot with name ' + title + " exists.")
        plot = self.addPlot(title = title)
        type = mesh.getROIType(roi_id)
        data = mesh.getROIData(roi_id)
        if type == steps.geom.ELEM_TRI:
            self.updater[title] = TrisSpecDistUpdater(plot, mesh, sim, data, spec_id, axis, nbins, y_range, **kwargs)
        elif type == steps.geom.ELEM_TET:
            self.updater[title] = TetsSpecDistUpdater(plot, mesh, sim, data, spec_id, axis, nbins, y_range, **kwargs)
        return plot
        
    def refresh(self):
        for p in self.updater.values():
            p.update()

    def reset(self):
        for p in self.updater.values():
            p.reset()
    

class CompSpecUpdater:

    def __init__(self, plot, sim, comp_id, spec_id, data_size, x_range = None, y_range = None, measure = "count", **kwargs):
        self.plot = plot
        plot.setLabel('bottom', 'Time', units='s')
        self.sim = sim
        self.comp_id = comp_id
        self.spec_id = spec_id
        self.measure = measure
        self.data_size = data_size
        self.time = [sim.getTime()]
        value = 0
        if self.measure == "count":
            value = sim.getCompCount(self.comp_id, self.spec_id)
        elif self.measure == "conc":
            value = sim.getCompConc(self.comp_id, self.spec_id)
        self.data = [value]
        self.curve = plot.plot(self.time, self.data, **kwargs)
        if x_range != None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])

    def update(self):
        self.time.append(self.sim.getTime())
        value = 0
        if self.measure == "count":
            value = self.sim.getCompCount(self.comp_id, self.spec_id)
        elif self.measure == "conc":
            value = self.sim.getCompConc(self.comp_id, self.spec_id)
        self.data.append(value)
        if len(self.time) >self.data_size:
            self.time.pop(0)
            self.data.pop(0)
        self.curve.setData(self.time, self.data)

    def reset(self):
        self.time = [self.sim.getTime()]
        value = 0
        if self.measure == "count":
            value = self.sim.getCompCount(self.comp_id, self.spec_id)
        elif self.measure == "conc":
            value = self.sim.getCompConc(self.comp_id, self.spec_id)
        self.data = [value]
        self.curve.setData(self.time, self.data)

class TetsSpecUpdater:
    
    def __init__(self, plot, sim, tets, spec_id, data_size, x_range = None, y_range = None, measure = "count", **kwargs):
        self.plot = plot
        plot.setLabel('bottom', 'Time', units='s')
        self.sim = sim
        self.tets = tets
        self.spec_id = spec_id
        self.measure = measure
        self.data_size = data_size
        self.time = [sim.getTime()]
        value = 0
        for t in self.tets:
            value += sim.getTetCount(t, self.spec_id)
        self.data = [value]
        self.curve = plot.plot(self.time, self.data, **kwargs)
        if x_range != None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
    
    def update(self):
        self.time.append(self.sim.getTime())
        value = 0
        for t in self.tets:
            value += self.sim.getTetCount(t, self.spec_id)
        self.data.append(value)
        if len(self.time) >self.data_size:
            self.time.pop(0)
            self.data.pop(0)
        self.curve.setData(self.time, self.data)
    
    def reset(self):
        self.time = [self.sim.getTime()]
        value = 0
        for t in self.tets:
            value += self.sim.getTetCount(t, self.spec_id)
        self.data = [value]
        self.curve.setData(self.time, self.data)

class PatchSpecUpdater:
    
    def __init__(self, plot, sim, patch_id, spec_id, data_size, x_range = None, y_range = None, **kwargs):
        self.plot = plot
        plot.setLabel('bottom', 'Time', units='s')
        self.sim = sim
        self.patch_id = patch_id
        self.spec_id = spec_id
        self.data_size = data_size
        self.time = [sim.getTime()]
        self.data = [sim.getPatchCount(patch_id, spec_id)]
        self.curve = plot.plot(self.time, self.data, **kwargs)
        if x_range != None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
    
    def update(self):
        self.time.append(self.sim.getTime())
        self.data.append(self.sim.getPatchCount(self.patch_id, self.spec_id))
        if len(self.time) >self.data_size:
            self.time.pop(0)
            self.data.pop(0)
        self.curve.setData(self.time, self.data)
    
    def reset(self):
        self.time = [self.sim.getTime()]
        self.data = [self.sim.getPatchCount(self.patch_id, self.spec_id)]
        self.curve.setData(self.time, self.data)

class TrisSpecUpdater:
    
    def __init__(self, plot, sim, tris, spec_id, data_size, x_range = None, y_range = None, **kwargs):
        self.plot = plot
        plot.setLabel('bottom', 'Time', units='s')
        self.sim = sim
        self.tris = tris
        self.spec_id = spec_id
        self.measure = measure
        self.data_size = data_size
        self.time = [sim.getTime()]
        value = 0
        for t in self.tris:
            value += sim.getTriCount(t, self.spec_id)
        self.data = [value]
        self.curve = plot.plot(self.time, self.data, **kwargs)
        if x_range != None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
    
    def update(self):
        self.time.append(self.sim.getTime())
        value = 0
        for t in self.tris:
            value += self.sim.getTriCount(t, self.spec_id)
        self.data.append(value)
        if len(self.time) >self.data_size:
            self.time.pop(0)
            self.data.pop(0)
        self.curve.setData(self.time, self.data)
    
    def reset(self):
        self.time = [self.sim.getTime()]
        value = 0
        for t in self.tris:
            value += self.sim.getTriCount(t, self.spec_id)
        self.data = [value]
        self.curve.setData(self.time, self.data)

class CompSumSpecsUpdater:
    
    def __init__(self, plot, sim, comp_id, spec_ids, data_size, x_range = None, y_range = None, **kwargs):
        self.plot = plot
        plot.setLabel('bottom', 'Time', units='s')
        self.sim = sim
        self.comp_id = comp_id
        self.spec_ids = spec_ids
        self.data_size = data_size
        self.time = [sim.getTime()]
        sum = 0
        for s in self.spec_ids:
            sum += sim.getCompCount(comp_id, s)
        self.data = [sum]
        self.curve = plot.plot(self.time, self.data, **kwargs)
        if x_range != None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
    
    def update(self):
        self.time.append(self.sim.getTime())
        sum = 0
        for s in self.spec_ids:
            sum += self.sim.getCompCount(self.comp_id, s)
        self.data.append(sum)
        if len(self.time) >self.data_size:
            self.time.pop(0)
            self.data.pop(0)
        self.curve.setData(self.time, self.data)
    
    def reset(self):
        self.time = [self.sim.getTime()]
        sum = 0
        for s in self.spec_ids:
            sum += self.sim.getCompCount(self.comp_id, s)
        self.data = [sum]
        self.curve.setData(self.time, self.data)

class PatchSumSpecsUpdater:
    
    def __init__(self, plot, sim, patch_id, spec_ids, data_size, x_range = None, y_range = None, **kwargs):
        self.plot = plot
        plot.setLabel('bottom', 'Time', units='s')
        self.sim = sim
        self.patch_id = patch_id
        self.spec_ids = spec_ids
        self.data_size = data_size
        self.time = [sim.getTime()]
        sum = 0
        for s in self.spec_ids:
            sum += sim.getPatchCount(patch_id, s)
        self.data = [sum]
        self.curve = plot.plot(self.time, self.data, **kwargs)
        if x_range != None:
            self.plot.setXRange(x_range[0], x_range[1])
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
    
    def update(self):
        self.time.append(self.sim.getTime())
        sum = 0
        for s in self.spec_ids:
            sum += self.sim.getPatchCount(self.patch_id, s)
        self.data.append(sum)
        if len(self.time) >self.data_size:
            self.time.pop(0)
            self.data.pop(0)
        self.curve.setData(self.time, self.data)
    
    def reset(self):
        self.time = [self.sim.getTime()]
        sum = 0
        for s in self.spec_ids:
            sum += self.sim.getPatchCount(self.patch_id, s)
        self.data = [sum]
        self.curve.setData(self.time, self.data)

class CompSpecDistUpdater:
    def __init__(self, plot, mesh, sim, comp_id, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        self.plot = plot
        self.mesh = mesh
        self.sim = sim
        self.comp_id = comp_id
        self.spec_id = spec_id
        self.nbins = nbins
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
        tmcomp = steps.geom.castToTmComp(mesh.getComp(comp_id))
        self.tets = tmcomp.getAllTetIndices()
        self.axis = 0
        if axis == "x":
            self.axis = 0
            plot.setLabel('bottom', 'X Coordinate', units='m')
        elif axis == "y":
            self.axis = 1
            plot.setLabel('bottom', 'Y Coordinate', units='m')
        elif axis == "z":
            self.axis = 2
            plot.setLabel('bottom', 'Z Coordinate', units='m')
        self.bound_max = tmcomp.getBoundMax()[self.axis]
        self.bound_min = tmcomp.getBoundMin()[self.axis]
        bins = np.linspace(self.bound_min, self.bound_max, nbins + 1)
        centers = [mesh.getTetBarycenter(t)[self.axis] for t in self.tets]
        self.belongs = np.digitize(centers, bins)
        self.bin_data = []
        for b in range(nbins):
            self.bin_data.append(bins[b])
            self.bin_data.append(bins[b + 1])
        data = np.zeros(nbins + 1)
        for t in range(len(self.tets)):
            data[self.belongs[t] - 1] += sim.getTetCount(self.tets[t], self.spec_id)
        y_data = []
        for b in range(nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        self.curve = plot.plot(self.bin_data, y_data, fillLevel = -0.3, **kwargs)


    def update(self):
        data = np.zeros(self.nbins + 1)
        for t in range(len(self.tets)):
            data[self.belongs[t] - 1] += self.sim.getTetCount(self.tets[t], self.spec_id)
        y_data = []
        for b in range(self.nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        
        self.curve.setData(self.bin_data, y_data)

    def reset(self):
        self.update()

class PatchSpecDistUpdater:
    def __init__(self, plot, mesh, sim, patch_id, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        self.plot = plot
        self.mesh = mesh
        self.sim = sim
        self.patch_id = patch_id
        self.spec_id = spec_id
        self.nbins = nbins
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
        tmpatch = steps.geom.castToTmPatch(mesh.getPatch(patch_id))
        self.tris = tmpatch.getAllTriIndices()
        self.axis = 0
        if axis == "x":
            self.axis = 0
            plot.setLabel('bottom', 'X Coordinate', units='m')
        elif axis == "y":
            self.axis = 1
            plot.setLabel('bottom', 'Y Coordinate', units='m')
        elif axis == "z":
            self.axis = 2
            plot.setLabel('bottom', 'Z Coordinate', units='m')
        self.bound_max = tmpatch.getBoundMax()[self.axis]
        self.bound_min = tmpatch.getBoundMin()[self.axis]
        bins = np.linspace(self.bound_min, self.bound_max, nbins + 1)
        centers = [mesh.getTriBarycenter(t)[self.axis] for t in self.tris]
        self.belongs = np.digitize(centers, bins)
        self.bin_data = []
        for b in range(nbins):
            self.bin_data.append(bins[b])
            self.bin_data.append(bins[b + 1])
        data = np.zeros(nbins + 1)
        for t in range(len(self.tris)):
            data[self.belongs[t] - 1] += sim.getTriCount(self.tris[t], self.spec_id)
        y_data = []
        for b in range(nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        self.curve = plot.plot(self.bin_data, y_data, fillLevel = -0.3, **kwargs)
    
    
    def update(self):
        data = np.zeros(self.nbins + 1)
        for t in range(len(self.tris)):
            data[self.belongs[t] - 1] += self.sim.getTriCount(self.tris[t], self.spec_id)
        y_data = []
        for b in range(self.nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        
        self.curve.setData(self.bin_data, y_data)
    
    def reset(self):
        self.update()
        
class TetsSpecDistUpdater:
    def __init__(self, plot, mesh, sim, tets, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        self.plot = plot
        self.mesh = mesh
        self.sim = sim
        self.spec_id = spec_id
        self.nbins = nbins
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
        self.tets = tets
        self.axis = 0
        if axis == "x":
            self.axis = 0
            plot.setLabel('bottom', 'X Coordinate', units='m')
        elif axis == "y":
            self.axis = 1
            plot.setLabel('bottom', 'Y Coordinate', units='m')
        elif axis == "z":
            self.axis = 2
            plot.setLabel('bottom', 'Z Coordinate', units='m')
        self.bound_max = -float("Inf")
        self.bound_min = float("Inf")
        centers = []
        for t in tets:
            verts = mesh.getTet(t)
            for v in verts:
                coordinate = mesh.getVertex(v)[self.axis]
                if coordinate > self.bound_max:
                    self.bound_max = coordinate
                if coordinate < self.bound_min:
                    self.bound_min = coordinate
            centers.append(mesh.getTetBarycenter(t)[self.axis])
        bins = np.linspace(self.bound_min, self.bound_max, nbins + 1)
        self.belongs = np.digitize(centers, bins)
        self.bin_data = []
        for b in range(nbins):
            self.bin_data.append(bins[b])
            self.bin_data.append(bins[b + 1])
        data = np.zeros(nbins + 1)
        for t in range(len(self.tets)):
            data[self.belongs[t] - 1] += sim.getTetCount(self.tets[t], self.spec_id)
        y_data = []
        for b in range(nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        self.curve = plot.plot(self.bin_data, y_data, fillLevel = -0.3, **kwargs)
    
    
    def update(self):
        data = np.zeros(self.nbins + 1)
        for t in range(len(self.tets)):
            data[self.belongs[t] - 1] += self.sim.getTetCount(self.tets[t], self.spec_id)
        y_data = []
        for b in range(self.nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        
        self.curve.setData(self.bin_data, y_data)
    
    def reset(self):
        self.update()

class TrisSpecDistUpdater:
    def __init__(self, plot, mesh, sim, tris, spec_id, axis = "x", nbins = 20, y_range = None, **kwargs):
        self.plot = plot
        self.mesh = mesh
        self.sim = sim
        self.spec_id = spec_id
        self.nbins = nbins
        if y_range != None:
            self.plot.setYRange(y_range[0], y_range[1])
        self.tris = tris
        self.axis = 0
        if axis == "x":
            self.axis = 0
            plot.setLabel('bottom', 'X Coordinate', units='m')
        elif axis == "y":
            self.axis = 1
            plot.setLabel('bottom', 'Y Coordinate', units='m')
        elif axis == "z":
            self.axis = 2
            plot.setLabel('bottom', 'Z Coordinate', units='m')
        self.bound_max = -float("Inf")
        self.bound_min = float("Inf")
        centers = []
        for t in tris:
            verts = mesh.getTri(t)
            for v in verts:
                coordinate = mesh.getVertex(v)[self.axis]
                if coordinate > self.bound_max:
                    self.bound_max = coordinate
                if coordinate < self.bound_min:
                    self.bound_min = coordinate
            centers.append(mesh.getTriBarycenter(t)[self.axis])
        bins = np.linspace(self.bound_min, self.bound_max, nbins + 1)
        self.belongs = np.digitize(centers, bins)
        self.bin_data = []
        for b in range(nbins):
            self.bin_data.append(bins[b])
            self.bin_data.append(bins[b + 1])
        data = np.zeros(nbins + 1)
        for t in range(len(self.tets)):
            data[self.belongs[t] - 1] += sim.getTriCount(self.tris[t], self.spec_id)
        y_data = []
        for b in range(nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        self.curve = plot.plot(self.bin_data, y_data, fillLevel = -0.3, **kwargs)
    
    
    def update(self):
        data = np.zeros(self.nbins + 1)
        for t in range(len(self.tris)):
            data[self.belongs[t] - 1] += self.sim.getTriCount(self.tris[t], self.spec_id)
        y_data = []
        for b in range(self.nbins):
            y_data.append(data[b])
            y_data.append(data[b])
        
        self.curve.setData(self.bin_data, y_data)
    
    def reset(self):
        self.update()



