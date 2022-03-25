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

from pyqtgraph.Qt import QtCore, QtGui, QtOpenGL
import pyqtgraph.opengl as gl
import numpy as np
import random

import steps.API_1.geom as sgeom
from steps.API_1.geom import INDEX_DTYPE
from steps.API_1.geom import UNKNOWN_TET
from steps.API_1.geom import UNKNOWN_TRI
import steps.API_1.utilities.meshctrl as smeshctrl

class VisualCompMesh(gl.GLMeshItem):
    """
    Static mesh component for a compartment
    Parameters:

    * id                      ID of the component
    * display                 Parent display
    * steps_mesh              STEPS mesh
    * comp_id                 ID of the compartment
    * color                   Color of the component
    """
    def __init__(self, id, display, steps_mesh, comp_id, color = None):
        """
        Constructor
        """
        self.id = id
        self. display = display
        self.steps_mesh = steps_mesh
        
        if not color:
            color = [random.random(), random.random(), random.random(), 0.3]
        tmcomp = sgeom.castToTmComp(steps_mesh.getComp(comp_id))
        self.bound_max = [v * display.scale for v in tmcomp.getBoundMax()]
        self.bound_min = [v * display.scale for v in tmcomp.getBoundMin()]

        ipatches = set(tmcomp.getIPatches())
        opatches = set(tmcomp.getOPatches())
        patches = ipatches.union(opatches)

        surface_tris = smeshctrl.findSurfTrisInComp(steps_mesh, tmcomp)
        for p in patches:
            surface_tris.extend(sgeom.castToTmPatch(p).getAllTriIndices())
        
        surface_tris = np.array(surface_tris, dtype = INDEX_DTYPE)
        
        
        v_set_size = steps_mesh.getTriVerticesSetSizeNP(surface_tris)
        tris_data = np.zeros(surface_tris.size * 3, dtype = INDEX_DTYPE)
        v_set = np.zeros(v_set_size, dtype = INDEX_DTYPE)
        verts_data = np.zeros(v_set_size * 3)
        steps_mesh.getTriVerticesMappingSetNP(surface_tris, tris_data, v_set)
        steps_mesh.getBatchVerticesNP(v_set, verts_data)
        verts_data *= display.scale
        tris_data.shape = -1, 3
        verts_data.shape = -1, 3
        mesh_data = gl.MeshData(vertexes=verts_data, faces = tris_data)
        gl.GLMeshItem.__init__(self, meshdata=mesh_data, smooth=False, computeNormals =True, shader='balloon', glOptions='additive')
        self.setColor(color)
        display.addItem(self)

    def updateItem(self):
        return

class VisualPatchMesh(gl.GLMeshItem):
    """
    Static mesh component for a patch
    
    Parameters:

    * id                      ID of the component
    * display                 Parent display
    * steps_mesh              STEPS mesh
    * patch_id                ID of the patch
    * color                   Color of the component
    """
    def __init__(self, id, display, steps_mesh, patch_id, color = None):
        """
        Constructor
        """
        self.id = id
        self. display = display
        self.steps_mesh = steps_mesh
        
        if not color:
            color = [random.random(), random.random(), random.random(), 0.3]
        
        tmpatch = sgeom.castToTmPatch(steps_mesh.getPatch(patch_id))
        self.bound_max = [v * display.scale for v in tmpatch.getBoundMax()]
        self.bound_min = [v * display.scale for v in tmpatch.getBoundMin()]
        
        patch_surface = tmpatch.getAllTriIndices()
        patch_surface = np.array(patch_surface, dtype = INDEX_DTYPE)
        v_set_size = steps_mesh.getTriVerticesSetSizeNP(patch_surface)
        tris_data = np.zeros(patch_surface.size * 3, dtype = INDEX_DTYPE)
        v_set = np.zeros(v_set_size, dtype = INDEX_DTYPE)
        verts_data = np.zeros(v_set_size * 3)
        steps_mesh.getTriVerticesMappingSetNP(patch_surface, tris_data, v_set)
        steps_mesh.getBatchVerticesNP(v_set, verts_data)
        verts_data *= display.scale
        tris_data.shape = -1, 3
        verts_data.shape = -1, 3
        mesh_data = gl.MeshData(vertexes=verts_data, faces = tris_data)
        gl.GLMeshItem.__init__(self, meshdata=mesh_data, smooth=False, computeNormals =True, shader='balloon', glOptions='additive')
        self.setColor(color)
        display.addItem(self)

    def updateItem(self):
        return



