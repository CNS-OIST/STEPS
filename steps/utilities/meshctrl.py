# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2011 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev: 352 $
#  Last Changed Date: $Date: 2010-08-09 10:49:54 +0900 (Mon, 09 Aug 2010) $
#  Last Changed By:   $Author: wchen $

""" 
.. Note::

    This module is preliminary, means some of the functions are still under development.
    Code modification / debugging is wellcomed.
    Please email steps.dev@gmail.com if you would like to share you changes with others.

Mesh Control Utilities 

The meshctrl module provides functions for controling and manipulating the mesh, 
e.g. finding overlap surface triangles in meshes.

"""

import steps.geom
from math import *

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################
def findOverlapTris(mesh, tets1, tets2):
    """
    Find overlap triangles between two sets of tetrahedrons within a mesh.
    
    Arguements:
        * steps.geom.Tetmesh mesh
        * list<uint> tets1
        * list<uint> tets2
        
    Return:
        list<uint>
    """
    
    tris1 = set()
    for i in tets1:
        tritemp = mesh.getTetTriNeighb(i)
        for j in range(4): 
            tris1.add(tritemp[j])
    tris2 = set()
    for i in tets2:
        tritemp = mesh.getTetTriNeighb(i)
        for j in range(4): 
            tris2.add(tritemp[j])

    common_tris = tris1.intersection(tris2)
    return list(common_tris)
    
################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################

def findOverlapSurfTris(mesh1, mesh2):
    """
    Find overlap surface triangles between two meshes.
    Return a list of coupling data list formatted as
    
    Arguements:
        * steps.geom.Tetmesh mesh1
        * steps.geom.Tetmesh mesh2
        
    Return:
        list<list<uint, uint, uint, uint, float>>
    """
    
    surfs1 = mesh1.getSurfTris()
    surfs2 = mesh2.getSurfTris()
    data = []
    
    for s1 in surfs1:
        c1 = mesh1.getTriBarycenter(s1)
        
        neighb_tets1 = mesh1.getTriTetNeighb(s1)
        tetid1 = neighb_tets1[0]
        if tetid1 < 0:
            tetid1 = neighb_tets1[1]
        assert(tetid1 >= 0)
        tc1 = mesh1.getTetBarycenter(tetid1)
        
        for s2 in surfs2:
            c2 = mesh2.getTriBarycenter(s2)
            if (c1[0] == c2[0]) and (c1[1] == c2[1]) and (c1[2] == c2[2]): 
                neighb_tets2 = mesh2.getTriTetNeighb(s2)
                tetid2 = neighb_tets2[0]
                if tetid2 < 0:
                    tetid2 = neighb_tets2[1]
                assert(tetid2 >= 0)
                tc2 = mesh2.getTetBarycenter(tetid2)
                
                xdist = tc1[0] - tc2[0]
                ydist = tc1[1] - tc2[1]
                zdist = tc1[2] - tc2[2]
                dist = sqrt((xdist*xdist) + (ydist*ydist) + (zdist*zdist))
                
                coupling = [tetid1, s1, tetid2, s2, dist]
                data.append(coupling)
    return data
                
################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################

def findSurfTrisInComp(mesh, comp_id):
    """
    Find surface triangles within a compartment.
    Return a list of surface triangle indices.
    
    Arguements:
        * steps.geom.Tetmesh mesh
        * string comp_id
        
    Return:
        list<uint>
    """
    
    surf_tris = mesh.getSurfTris()
    comp = mesh.getComp(comp_id)
    comp_tets = comp.getAllTetIndices()
    
    out_surfs = []
    for tet in  comp_tets:
        nb_tris = mesh.getTriTetNeighb(tet)
        for tri in nb_tris:
            if tri in surf_tris:
                out_surfs.append(tri)
    
    return out_surfs
            
################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################

def findSurfTrisInTets(mesh, tet_list):
    """
    Find surface triangles within a compartment.
    Return a list of surface triangle indices.
    
    Arguements:
        * steps.geom.Tetmesh mesh
        * string comp_id
        
    Return:
        list<uint>
    """
    
    surf_tris = mesh.getSurfTris()

    out_surfs = []
    for tet in  tet_list:
        nb_tris = mesh.getTriTetNeighb(tet)
        for tri in nb_tris:
            if tri in surf_tris:
                out_surfs.append(tri)
    
    return out_surfs
            
################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################