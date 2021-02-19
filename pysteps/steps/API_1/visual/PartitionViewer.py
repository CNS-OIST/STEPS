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
from __future__ import print_function

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtOpenGL
import numpy as np
import pyqtgraph.opengl as gl
import random
from steps.API_1.geom import INDEX_DTYPE
from steps.API_1.geom import UNKNOWN_TET
from steps.API_1.geom import UNKNOWN_TRI
def createColorMap(partitions):
    """
    Genrate a random color map for a partition.
    Parameter:

    * partitions           Partitions for a list of elements. The partition can be either a list, or a dict.
       
    Return:
        Return a dict color_map, for each partion in partitions, color_map[partition] gives a randomly generated color.
    """

    if isinstance(partitions, list):
        color_map = {}
        for part in partitions:
            if part not in color_map.keys():
                if part is None:
                    color_map[part] = [1.0, 0.0, 0.0, 1.0]
                else:
                    color_map[part] = [random.random(), random.random(), random.random(), 0.3]
        return color_map
    elif isinstance(partitions, dict):
        color_map = {}
        for part in partitions.values():
            if part not in color_map.keys():
                if part is None:
                    color_map[part] = [1.0, 0.0, 0.0, 1.0]
                else:
                    color_map[part] = [random.random(), random.random(), random.random(), 0.3]
        return color_map

    print("Unknown partitioning data, no color map is generated.")
    return None

class TetPartitionDisplay(QtGui.QMainWindow):
    """
    Partition Display
    Parameters:

    * mesh                STEPS Tetmesh
    * tet_partitions      Partition of tetrahedrons
    * title               Display title
    * x                   X cordinate of the display
    * y                   Y cordinate of the display
    * w                   Width of the display
    * h                   Height of the display
    * scale               Scaling between STEPS mesurement and display pixel
    * color_map           Color map for each partition, it is a Python dict with key-value pairs as color_map[partition_id] = [red, green, blue, alpha], where partition_id is the data stored in tet_partitions. The color is defined by four parameters [red, green, blue, alpha], each with range from 0.0 to 1.0. If the partition name is not in the color map a random color will be generated. A partition assigned to None will always be colored as nontransparent red [1.0, 0.0, 0.0, 1.0].
    * morph_sections        NEURON hoc morphology data of the mesh
    * show_list           List of sections to be visualized
    """
    
    def __init__(self, mesh, tet_partitions, title = "TetPartitionDisplay", x = 100, y = 100, w = 800, h = 600, scale = 1e6, color_map = None, morph_sections = None,  show_list = None):
        """
        Constructor.
        """
        QtGui.QMainWindow.__init__(self)
        self.setGeometry(x, y, w, h)
        self.widget = gl.GLViewWidget(self)
        self.setCentralWidget(self.widget)
        self.widget_mapping = {}
        self.mesh = mesh
        self.tet_part_table = {}
        for tet in range(len(tet_partitions)):
            if tet_partitions[tet] not in self.tet_part_table.keys():
                self.tet_part_table[tet_partitions[tet]] = []
            self.tet_part_table[tet_partitions[tet]].append(tet)
        
        self.scale = scale
        self.setWindowTitle(title)
        
        self.bound_min = [v * self.scale for v in mesh.getBoundMin()]
        self.bound_max = [v * self.scale for v in mesh.getBoundMax()]
        
        self.center = [0.0, 0.0, 0.0]
        
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
        if color_map is None:
            self.color_map = {}
        else:
            self.color_map = color_map
        for tet_part in self.tet_part_table:
            if tet_part is None:
                self.color_map[None] = [1.0, 0.0, 0.0, 1.0]
            else:
                if tet_part not in self.color_map:
                    color = [random.random(), random.random(), random.random(), 0.3]
                    self.color_map[tet_part] = color
            part = TetPartitionMesh(self, mesh, self.tet_part_table[tet_part], color = self.color_map[tet_part])
            self.widget.addItem(part)
            if tet_part not in self.widget_mapping:
                self.widget_mapping[tet_part] = []
            self.widget_mapping[tet_part].append(part)

        if morph_sections is not None:
            for sec in morph_sections.values():
                c = None
                if sec["name"] not in self.color_map:
                    print("Section ", sec["name"], " does not have defined color, generate random one instead.")
                    c = [random.random(), random.random(), random.random(), 0.3]
                    self.color_map[sec["name"]] = c
                else:
                    c = self.color_map[sec["name"]]
                data = []
                points = sec["points"]
                npoints = len(points)
                for i in range(npoints - 1):
                    p0 = points[i]
                    p1 = points[i+1]
                    data.append([p0[0], p0[1], p0[2]])
                    data.append([p1[0], p1[1], p1[2]])
                data = np.array(data)
                data.shape = -1, 3

                line_color = [c[0], c[1], c[2], 1.0]
                section = gl.GLLinePlotItem(pos = data, color = line_color, width = 3, mode = "lines")
                self.widget.addItem(section)
                if sec["name"] not in self.widget_mapping:
                    self.widget_mapping[sec["name"]] = []
                self.widget_mapping[sec["name"]].append(section)
        self.show()

    def showSections(self, section_list):
        for w in self.widget.items:
            w.hide()
        for sec in section_list:
            if sec not in self.widget_mapping:
                print("Section ", sec, " does not have visualization compoment.")
                continue
            for w in self.widget_mapping[sec]:
                w.show()

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

class TetPartitionMesh(gl.GLMeshItem):
    """
    Static mesh component for a compartment
    Parameters:

    * display                 Parent display
    * steps_mesh              STEPS mesh
    * tet_list                Tetrahedron list of a section of the mesh
    * color                   Color of the component
    """
    def __init__(self, display, mesh, tet_list, color = None):
        """
        Constructor
        """
        self.display = display
        self.mesh = mesh
        
        if not color:
            color = [random.random(), random.random(), random.random(), 0.3]
        
        surface_tris = []
        
        for tet in tet_list:
            tris = mesh.getTetTriNeighb(tet)
            for tri in tris:
                if tri == UNKNOWN_TRI: continue
                neighb_tets = mesh.getTriTetNeighb(tri)
                for neighb_tet in neighb_tets:
                    if neighb_tet == UNKNOWN_TET or neighb_tet not in tet_list:
                        surface_tris.append(tri)
    
        surface_tris = np.array(surface_tris, dtype = INDEX_DTYPE)
        
        v_set_size = mesh.getTriVerticesSetSizeNP(surface_tris)
        tris_data = np.zeros(surface_tris.size * 3, dtype = INDEX_DTYPE)
        v_set = np.zeros(v_set_size, dtype = INDEX_DTYPE)
        verts_data = np.zeros(v_set_size * 3)
        mesh.getTriVerticesMappingSetNP(surface_tris, tris_data, v_set)
        mesh.getBatchVerticesNP(v_set, verts_data)
        verts_data *= display.scale
        tris_data.shape = -1, 3
        verts_data.shape = -1, 3
        mesh_data = gl.MeshData(vertexes=verts_data, faces = tris_data)
        gl.GLMeshItem.__init__(self, meshdata=mesh_data, smooth=False, computeNormals =True, shader='balloon', glOptions='additive')
        self.setColor(color)

class TriPartitionDisplay(QtGui.QMainWindow):
    """
    Partition Display
    Parameters:

    * mesh                STEPS Tetmesh
    * tri_partitions      Partition of triangles
    * title               Display title
    * x                   X cordinate of the display
    * y                   Y cordinate of the display
    * w                   Width of the display
    * h                   Height of the display
    * scale               Scaling between STEPS mesurement and display pixel
    * color_map           Color map for each partition, it is a Python dict with key-value pairs as color_map[partition_id] = [red, green, blue, alpha], where partition_id is the data stored in tet_partitions. The color is defined by four parameters [red, green, blue, alpha], each with range from 0.0 to 1.0. If the partition name is not in the color map a random color will be generated. A partition assigned to None will always be colored as nontransparent red [1.0, 0.0, 0.0, 1.0].
    * morph_sections      NEURON hoc morphology data of the mesh
    * show_list           List of sections to be visualized
    """
    
    def __init__(self, mesh, tri_partitions, title = "TriPartitionDisplay", x = 100, y = 100, w = 800, h = 600, scale = 1e6, color_map = None, morph_sections = None, show_list = None):
        """
        Constructor.
        """
        QtGui.QMainWindow.__init__(self)
        self.setGeometry(x, y, w, h)
        self.widget = gl.GLViewWidget(self)
        self.setCentralWidget(self.widget)
        self.widget_mapping = {}
        self.mesh = mesh
        
        self.tri_part_table = {}
        for tri in tri_partitions:
            if tri_partitions[tri] not in self.tri_part_table.keys():
                self.tri_part_table[tri_partitions[tri]] = []
            self.tri_part_table[tri_partitions[tri]].append(tri)
        
        self.scale = scale
        self.setWindowTitle(title)
        
        self.bound_min = [v * self.scale for v in mesh.getBoundMin()]
        self.bound_max = [v * self.scale for v in mesh.getBoundMax()]
        
        self.center = [0.0, 0.0, 0.0]
        
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
        if color_map == None:
            self.color_map = {}
        else:
            self.color_map = color_map
        for tri_part in self.tri_part_table:
            if tri_part == None:
                self.color_map[None] = [1.0, 0.0, 0.0, 1.0]
            else:
                if tri_part not in self.color_map:
                    color = [random.random(), random.random(), random.random(), 0.3]
                    self.color_map[tri_part] = color
            part = TriPartitionMesh(self, mesh, self.tri_part_table[tri_part], color = self.color_map[tri_part])
            self.widget.addItem(part)
            if tri_part not in self.widget_mapping:
                self.widget_mapping[tri_part] = []
            self.widget_mapping[tri_part].append(part)

        if morph_sections != None:
            for sec in morph_sections.values():
                c = None
                if sec["name"] not in self.color_map:
                    print("Section ", sec["name"], " does not have defined color, generate random one instead.")
                    c = [random.random(), random.random(), random.random(), 0.3]
                    self.color_map[sec["name"]] = c
                else:
                    c = self.color_map[sec["name"]]
                data = []
                points = sec["points"]
                npoints = len(points)
                for i in range(npoints - 1):
                    p0 = points[i]
                    p1 = points[i+1]
                    data.append([p0[0], p0[1], p0[2]])
                    data.append([p1[0], p1[1], p1[2]])
                data = np.array(data)
                data.shape = -1, 3

                line_color = [c[0], c[1], c[2], 1.0]
                section = gl.GLLinePlotItem(pos = data, color = line_color, width = 3, mode = "lines")
                self.widget.addItem(section)
                if sec["name"] not in self.widget_mapping:
                    self.widget_mapping[sec["name"]] = []
                self.widget_mapping[sec["name"]].append(section)
        self.show()

    def showSections(self, section_list):
        for w in self.widget.items:
            w.hide()
        for sec in section_list:
            if sec not in self.widget_mapping:
                print("Section ", sec, " does not have visualization compoment.")
                continue
            for w in self.widget_mapping[sec]:
                w.show()


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

class TriPartitionMesh(gl.GLMeshItem):
    """
    Static mesh component for a compartment
    Parameters:

    * display                 Parent display
    * steps_mesh              STEPS mesh
    * tri_list                Triangle list of a section of the mesh
    * color                   Color of the component
    """
    def __init__(self, display, mesh, tri_list, color = None):
        """
        Constructor
        """
        self.display = display
        self.mesh = mesh
        
        if not color:
            color = [random.random(), random.random(), random.random(), 0.3]
    
        surface_tris = np.array(tri_list, dtype = INDEX_DTYPE)
        
        v_set_size = mesh.getTriVerticesSetSizeNP(surface_tris)
        tris_data = np.zeros(surface_tris.size * 3, dtype = INDEX_DTYPE)
        v_set = np.zeros(v_set_size, dtype = INDEX_DTYPE)
        verts_data = np.zeros(v_set_size * 3)
        mesh.getTriVerticesMappingSetNP(surface_tris, tris_data, v_set)
        mesh.getBatchVerticesNP(v_set, verts_data)
        verts_data *= display.scale
        tris_data.shape = -1, 3
        verts_data.shape = -1, 3
        mesh_data = gl.MeshData(vertexes=verts_data, faces = tris_data)
        gl.GLMeshItem.__init__(self, meshdata=mesh_data, smooth=False, computeNormals =True, shader='balloon', glOptions='additive')
        self.setColor(color)


