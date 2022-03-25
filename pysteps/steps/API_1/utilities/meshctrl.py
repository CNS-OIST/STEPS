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

""" 
.. Note::

    This module is preliminary, means some of the functions are still under development.
    Code modification / debugging is wellcomed.
    Please email steps.dev@gmail.com if you would like to share you changes with others.

Mesh Control Utilities 

The meshctrl module provides functions for controling and manipulating the mesh, 
e.g. finding overlap surface triangles in meshes.

"""

import steps.API_1.geom
from steps.API_1.geom import UNKNOWN_TET
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
    [tet_1, surftri_1, tet_2, surftri_2, tet_distance],
    where:
    * tet_1 and tet_2: Indices of the tetrahedrons 
    whose surface triangles overlap.
    * surftri_1 and surftri_2: Indices of the overlap surface triangles.
    * tet_distance: Distance between the barycenters of tet_1 and tet_2.
    
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
        if tetid1 == UNKNOWN_TET:
            tetid1 = neighb_tets1[1]
        assert(tetid1 >= 0)
        tc1 = mesh1.getTetBarycenter(tetid1)
        
        for s2 in surfs2:
            c2 = mesh2.getTriBarycenter(s2)
            if (c1[0] == c2[0]) and (c1[1] == c2[1]) and (c1[2] == c2[2]): 
                neighb_tets2 = mesh2.getTriTetNeighb(s2)
                tetid2 = neighb_tets2[0]
                if tetid2 == UNKNOWN_TET:
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

def findSurfTrisInComp(mesh, comp):
    """
        Find surface triangles within a compartment.
        Return a list of surface triangle indices.
        
        Arguements:
        * steps.geom.Tetmesh mesh
        * steps.geom.TmComp comp
        
        Return:
        list<uint>
    """
    
    comp_tets = comp.getAllTetIndices()

    return findSurfTrisInTets(mesh, comp_tets)
            
################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################

def findSurfTrisInTets(mesh, tet_list):
    """
    Find surface triangles within a list of tetrahedrons.
    Return a list of surface triangle indices.
    
    Arguements:
        * steps.geom.Tetmesh mesh
        * list<uint> tet_list
        
    Return:
        list<uint>
    """

    tet_tris = set()

    for tet in  tet_list:
        nb_tris = mesh.getTetTriNeighb(tet)
        for tri in nb_tris:
            tet_tris.add(tri)

    surf_tris = set(mesh.getSurfTris())

    return list(surf_tris.intersection(tet_tris))
            
################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################

