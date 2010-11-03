# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
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

try:
    import cubit
except ImportError:
    print("Unable to find CUBIT module.")
    
def getTetsInVolume(v_id, mesh, scale):
    ntets = mesh.ntets
    volume = cubit.volume(v_id)
    body = volume.bodies()[0]
    in_list = []
    for t in range(ntets):
        center = mesh.getTetBarycenter(t)
        cubit_center = [center[0]/scale, center[1]/scale, center[2]/scale]
        status = body.point_containment(cubit_center)
        if status == 1 or status == 2:
            in_list.append(t)
    return in_list
    
def getCompTetsInVolume(v_id, mesh, comp_id, scale):
    comp_tets = mesh.getComp(comp_id).getAllTetIndices()
    volume = cubit.volume(v_id)
    body = volume.bodies()[0]
    in_list = []
    for t in comp_tets:
        center = mesh.getTetBarycenter(t)
        cubit_center = [center[0]/scale, center[1]/scale, center[2]/scale]
        status = body.point_containment(cubit_center)
        if status == 1 or status == 2:
            in_list.append(t)
    return in_list
        
def getSurfTrisInVolume(v_id, mesh, scale):
    surf_tris = mesh.getSurfTris()
    volume = cubit.volume(v_id)
    body = volume.bodies()[0]
    in_list = []
    for t in surf_tris:
        center = mesh.getTriBarycenter(t)
        cubit_center = [center[0]/scale, center[1]/scale, center[2]/scale]
        status = body.point_containment(cubit_center)
        if status == 1 or status == 2:
            in_list.append(t)
    return in_list
    
def getPatchTrisInVolume(v_id, mesh, patch_id, scale):
    patch_tris = mesh.getPatch(patch_id).getAllTriIndices()
    volume = cubit.volume(v_id)
    body = volume.bodies()[0]
    in_list = []
    for t in patch_tris:
        center = mesh.getTriBarycenter(t)
        cubit_center = [center[0]/scale, center[1]/scale, center[2]/scale]
        status = body.point_containment(cubit_center)
        if status == 1 or status == 2:
            in_list.append(t)
    return in_list
    
def highlightTets(steps_tets, tet_proxy):
    cubit.cmd("Graphics Clear Highlight")
    cubit.cmd("Graphics Pause")
    cmd_str = "highlight tet "
    for st in steps_tets:
        cubit_id = tet_proxy.getImportID(st)
        cmd_str += "%i," % (cubit_id)
    cubit.cmd(cmd_str)
    cubit.cmd("Display")