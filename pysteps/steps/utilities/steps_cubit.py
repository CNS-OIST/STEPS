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
from __future__ import print_function

try:
    import cubit
except ImportError:
    print("Unable to import CUBIT module.")

from steps.utilities.steps_shadow import *

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getSelectedVolumes():
    """
    Return the CUBIT indices of the selected volumes.
        
    Parameters:
        None
        
    Return:
        cubit_ids
    """
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    idxs = cubit.get_group_volumes(group_id)
    cubit.delete_group(group_id)
    return idxs

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getSelectedSurfaces():
    """
    Return the CUBIT indices of the selected surfaces.
        
    Parameters:
        None
        
    Return:
        cubit_ids
    """
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    idxs = cubit.get_group_surfaces(group_id)
    cubit.delete_group(group_id)
    return idxs

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getSelectedNodes():
    """
    Return the CUBIT indices of the selected vertex nodes.
        
    Parameters:
        None
        
    Return:
        cubit_ids
    """
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    idxs = cubit.get_group_nodes(group_id)
    cubit.delete_group(group_id)
    return idxs

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getSelectedTets():
    """
    Return the CUBIT indices of the selected tetrahedrons.
    
    Parameters:
        None
        
    Return:
        cubit_ids
    """
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    idxs = cubit.get_group_tets(group_id)
    cubit.delete_group(group_id)
    return idxs

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getSelectedTris():
    """
    Return the CUBIT indices of the selected triangles.
    
    Parameters:
        None
    
    Return:
        cubit_ids
    """
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    idxs = cubit.get_group_tris(group_id)
    cubit.delete_group(group_id)
    return idxs

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedVolumesAsComp(comp_name, mesh, vsys):
    """
    Create a shadow compartment using selected volumes.
        
    Parameters:
        * comp_name   Name of the compartment
        * mesh        ShadowMesh object
        * vsys        List of STEPS volume system ids for the compartment
        
    Return:
        ShadowComp object
    """
    vols = getSelectedVolumes()
    idxs = []
    for v in vols:
        idxs.extend(cubit.get_volume_tets(v))
    return ShadowComp(comp_name, mesh, idxs, vsys)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedSurfacesAsPatch(patch_name, mesh, ssys, icomp, ocomp = None):
    """
    Create a shadow patch using selected surfaces.
        
    Parameters:
        * patch_name  Name of the patch
        * mesh        ShadowMesh object
        * ssys        List of STEPS volume system ids for the compartment
        * icomp       ShadowComp object as inner compartment
        * ocomp       ShadowComp object as outer compartment, None by default
        
    Return:
        ShadowPatch object
    """
    surfs = getSelectedSurfaces()
    idxs = []
    for s in surfs:
        idxs.extend(cubit.get_surface_tris(s))
    return ShadowPatch(patch_name, mesh, idxs, ssys, icomp, ocomp)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedTetsAsComp(comp_name, mesh, vsys):
    """
    Create a shadow compartment using selected tetraedrons.
    
    Parameters:
        * comp_name   Name of the compartment
        * mesh        ShadowMesh object
        * vsys        List of STEPS volume system ids for the compartment
        
    Return:
        ShadowComp object
    """
    
    idxs = getSelectedTets()
    return ShadowComp(comp_name, mesh, idxs, vsys)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedTrisAsPatch(patch_name, mesh, ssys, icomp, ocomp = None):
    """
    Create a shadow patch using selected triangles.
        
    Parameters:
        * patch_name  Name of the patch
        * mesh        ShadowMesh object
        * ssys        List of STEPS volume system ids for the compartment
        * icomp       ShadowComp object as inner compartment
        * ocomp       ShadowComp object as outer compartment, None by default
        
    Return:
        ShadowPatch object
    """
    
    idxs = getSelectedTris()
    return ShadowPatch(patch_name, mesh, idxs, ssys, icomp, ocomp)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedNodesAsROI(roi_name, mesh):
    """
    Store the selected nodes as Region of Interest in the ShadowMesh object.
        
    Parameters:
        * roi_name    Name of the ROI
        * mesh        ShadowMesh object

    Return:
        None
    """
    
    idxs = getSelectedNodes()
    mesh.addROI(roi_name, ELEM_VERTEX, idxs)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedTetsAsROI(roi_name, mesh):
    """
    Store the selected terrahedrons as Region of Interest in the ShadowMesh object.
        
    Parameters:
        * roi_name    Name of the ROI
        * mesh        ShadowMesh object

        
    Return:
        None
    """
    
    idxs = getSelectedTets()
    mesh.addROI(roi_name, ELEM_TET, idxs)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def selectedTrisAsROI(roi_name, mesh):
    """
    Store the selected triangles as Region of Interest in the ShadowMesh object.
        
    Parameters:
        * roi_name    Name of the ROI
        * mesh        ShadowMesh object

    Return:
        None
    """
    
    idxs = getSelectedTris()
    mesh.addROI(roi_name, ELEM_TRI, idxs)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getNodesBoundInSelectedVols(target_list):
    """
    Return nodes bound in selected CUBIT volumes.
        
    Parameters:
        * target_list List of indices of target nodes
                    
    Return:
        List of indices of nodes bound by the volumes
    """
    
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    volume_ids = cubit.get_group_volumes(group_id)
    
    in_list = []
    
    for v in target_list:
        cords = cubit.get_nodal_coordinates(v)
        for vol_id in volume_ids:
            volume = cubit.volume(vol_id)
            body = volume.bodies()[0]
            status = body.point_containment(cords)
            if status == 1 or status == 2:
                in_list.append(v)
                break
    cubit.delete_group(group_id)
    return in_list

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getTetsBoundInSelectedVols(target_list):
    """
    Return tetraherons bound in selected CUBIT volumes.
        
    Parameters:
        * target_list List of indices of target tetrahedrons
                    
    Return:
        List of indices of tetrahedrons bound by the volumes
    """
    
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    volume_ids = cubit.get_group_volumes(group_id)
    
    in_list = []
    
    for t in target_list:
        verts = cubit.get_connectivity("tet", t)
        with_in = True
        cords = []
        for v in verts:
            c = cubit.get_nodal_coordinates(v)
            cords.append(c)
        for vol_id in volume_ids:
            volume = cubit.volume(vol_id)
            body = volume.bodies()[0]
            within = True
            for cord in cords:
                status = body.point_containment(cord)
                if status == 0:
                    within = False
                    break
            if within:
                in_list.append(t)
                break
    cubit.delete_group(group_id)
    return in_list

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def getTrisBoundInSelectedVols(target_list):
    """
    Return triangles bound in selected CUBIT volumes.
        
    Parameters:
        * target_list List of indices of target triangles
                    
    Return:
        List of indices of triangles bound by the volumes
    """
    
    group_id = cubit.create_new_group()
    cubit.silent_cmd("group %i add selection" % (group_id))
    volume_ids = cubit.get_group_volumes(group_id)
    
    in_list = []
    
    for t in target_list:
        verts = cubit.get_connectivity("tri", t)
        with_in = True
        cords = []
        for v in verts:
            c = cubit.get_nodal_coordinates(v)
            cords.append(c)
        for vol_id in volume_ids:
            volume = cubit.volume(vol_id)
            body = volume.bodies()[0]
            within = True
            for cord in cords:
                status = body.point_containment(cord)
                if status == 0:
                    within = False
                    break
            if within:
                in_list.append(t)
                break
    cubit.delete_group(group_id)
    return in_list

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def boundTetsAsComp(target_list, comp_name, mesh, vsys):
    """
    Create a shadow compartment using tetraedrons bound in selected volumes.
        
    Parameters:
        * target_list List of indices of target tetrahedrons
        * comp_name   Name of the compartment
        * mesh        ShadowMesh object
        * vsys        List of STEPS volume system ids for the compartment
        
    Return:
        ShadowComp object
    """
    
    idxs = getTetsBoundInSelectedVols(target_list)
    return ShadowComp(comp_name, mesh, idxs, vsys)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def boundTrisAsPatch(target_list, patch_name, mesh, ssys, icomp, ocomp = None):
    """
    Create a shadow patch using triangles bound in selected volumes.
        
    Parameters:
        * target_list List of indices of target triangles
        * patch_name  Name of the patch
        * mesh        ShadowMesh object
        * ssys        List of STEPS volume system ids for the compartment
        * icomp       ShadowComp object as inner compartment
        * vocomp       ShadowComp object as outer compartment, None by default
        
    Return:
        ShadowPatch object
    """
    
    idxs = getTrisBoundInSelectedVols(target_list)
    return ShadowPatch(patch_name, mesh, idsx, ssys, icomp, ocomp)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def boundNodesAsROI(target_list, roi_name, mesh):
    """
    Store nodes bound in selected volumes as Region of Interest in the shadow mesh.
        
    Parameters:
        * target_list List of indices of target nodes
        * roi_name    Name of the ROI
        * mesh        ShadowMesh object
        
    Return:
        None
    """
    
    idxs = getNodesBoundInSelectedVols(target_list)
    mesh.addROI(roi_name, ELEM_VERTEX, idxs)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def boundTetsAsROI(target_list, roi_name, mesh):
    """
    Store tetrahedrons bound in selected volumes as Region of Interest in the shadow mesh.
        
    Parameters:
        * target_list List of indices of target tetrahedrons
        * roi_name    Name of the ROI
        * mesh        ShadowMesh object
        
    Return:
        None
    """
    
    idxs = getTetsBoundInSelectedVols(target_list)
    mesh.addROI(roi_name, ELEM_TET, idxs)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def boundTrisAsROI(target_list, roi_name, mesh):
    """
    Store triangles bound in selected volumes as Region of Interest in the shadow mesh.
        
    Parameters:
        * target_list List of indices of target triangles
        * roi_name    Name of the ROI
        * mesh        ShadowMesh object
        
    Return:
        None
    """
    
    idxs = getTrisBoundInSelectedVols(target_list)
    mesh.addROI(roi_name, ELEM_TRI, idxs)

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def drawROI(mesh, roi_name):
    """
    Draw elements stored in ROI database with id roi_name.
    
    Parameters:
        * mesh        ShadowMesh object
        * roi_name    Name of the ROI
        
    Return:
        None
    """
    if roi_name not in mesh.rois:
        raise Exception(roi_name + " is not in ROI database.")
    
    roi = mesh.rois[roi_name]
    idx_str = toStr(roi["Indices"])
    
    if roi["Type"] == ELEM_VERTEX:
        cubit.silent_cmd("draw node " + idx_str)
    elif roi["Type"] == ELEM_TET:
        cubit.silent_cmd("draw tet " + idx_str)
    elif roi["Type"] == ELEM_TRI:
        cubit.silent_cmd("draw tri " + idx_str)
    else:
        raise Exception(roi_name + " is undefined element type.")

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def drawComp(comp):
    """
    Draw the compartment.
    
    Parameters:
        * comp        ShadowComp object
        
    Return:
        None
    """
    cubit.silent_cmd("draw tet " + toStr(comp.indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def drawPatch(patch):
    """
    Draw the patch.
        
    Parameters:
        * patch       ShadowPatch object
        
    Return:
        None
    """
    cubit.silent_cmd("draw tri " + toStr(patch.indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def drawComp(mesh, comp_id):
    """
    Draw the compartment.
        
    Parameters:
        * mesh        ShadowMesh object
        * comp_id     ID of the ShadowComp object
        
    Return:
        None
    """
    cubit.silent_cmd("draw tet " + toStr(mesh.comps[comp_id].indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def drawPatch(mesh, patch_id):
    """
    Draw the patch.
        
    Parameters:
        * mesh        ShadowMesh object
        * patch_id    ID of the ShadowPatch object
        
    Return:
        None
    """
    cubit.silent_cmd("draw tri " + toStr(mesh.patches[patch_id].indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def highlightROI(mesh, roi_name):
    """
    Highlight elements stored in ROI database with id roi_name.
    
    Parameters:
        * mesh        ShadowMesh object
        * roi_name    Name of the ROI
        
    Return:
        None
    """
    if roi_name not in mesh.rois:
        raise Exception(roi_name + " is not in ROI database.")
    
    roi = mesh.rois[roi_name]
    idx_str = toStr(roi["Indices"])
    
    if roi["Type"] == ELEM_VERTEX:
        cubit.silent_cmd("highlight node " + idx_str)
    elif roi["Type"] == ELEM_TET:
        cubit.silent_cmd("highlight tet " + idx_str)
    elif roi["Type"] == ELEM_TRI:
        cubit.silent_cmd("highlight tri " + idx_str)
    else:
        raise Exception(roi_name + " is undefined element type.")

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def highlightComp(comp):
    """
    Highlight the compartment.
    
    Parameters:
        * comp        ShadowComp object
        
    Return:
        None
    """
    cubit.silent_cmd("highlight tet " + toStr(comp.indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def highlightPatch(patch):
    """
    Highlight the patch.
        
    Parameters:
        * patch       ShadowPatch object
        
    Return:
        None
    """
    cubit.silent_cmd("highlight tri " + toStr(patch.indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def highlightComp(mesh, comp_id):
    """
    Highlight the compartment.
        
    Parameters:
        * mesh        ShadowMesh object
        * comp_id     ID of the ShadowComp object
        
    Return:
        None
    """
    cubit.silent_cmd("highlight tet " + toStr(mesh.comps[comp_id].indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################

def highlightPatch(mesh, patch_id):
    """
    Highlight the patch.
        
    Parameters:
        * mesh        ShadowMesh object
        * patch_id    ID of the ShadowPatch object
        
    Return:
        None
    """
    cubit.silent_cmd("highlight tri " + toStr(mesh.patches[patch_id].indices))

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################

def toStr(e_list):
    """
    Convert entities list to a text string.
        
    Parameters:
        * e_list      entity index list
        
    Return:
        String of the entities separated by comma
    """
    return_str = ""
    for e in e_list:
        return_str += "%i," % (e)
    return return_str

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
