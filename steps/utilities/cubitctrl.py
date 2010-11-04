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
    
import random
import os
import sys

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################

def getTetsInVolume(v_id, mesh, scale = 1e-6):
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
 
################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################

def getCompTetsInVolume(v_id, mesh, comp_id, scale = 1e-6):
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

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################
        
def getSurfTrisInVolume(v_id, mesh, scale = 1e-6):
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
 
################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################

def getPatchTrisInVolume(v_id, mesh, patch_id, scale = 1e-6):
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

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################
    
def highlightTets(steps_tets, tet_proxy):
    cubit.cmd("Graphics Clear Highlight")
    cubit.cmd("Graphics Pause")
    cmd_str = "highlight tet "
    for st in steps_tets:
        cubit_id = tet_proxy.getImportID(st)
        cmd_str += "%i," % (cubit_id)
    cubit.cmd(cmd_str)
    cubit.cmd("Display")

################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################
    
class VSim(object):
    def __init__(self, model, mesh, sim, scale = 1e-6,color_map = None,
                cubit_file = None, visual_mode = "wireframe", with_mesh = False):
        if cubit_file != None:
            cubit.cmd('open "%s"' % (cubit_file))
        
        cubit.cmd(visual_mode)
        if with_mesh:
            cubit.cmd("mesh on")
        else:
            cubit.cmd("mesh off")
        self.model = model
        self.mesh = mesh
        self.sim = sim
        self.scale = scale

        self.specs = []
        spec_refs = model.getAllSpecs()
        for spec in spec_refs:
            self.specs.append(spec.getID())

        colors = random.sample(range(86), len(self.specs))

        self.colors = {}
        for i in range(len(colors)):
            self.colors[self.specs[i]] = colors[i]
            
        if color_map != None:
            for s in color_map:
                self.colors[s] = color_map[s]
                
        comp_refs = mesh.getAllComps()
        self.comp_tets = {}
        for c in comp_refs:
            self.comp_tets[c.getID()] = []
        for t in range(mesh.ntets):
            comp = self.mesh.getTetComp(t)
            if comp != None:
                self.comp_tets[comp.getID()].append(t)
                
        patch_refs = mesh.getAllPatches()
        self.patch_tris = {}
        for p in patch_refs:
            self.patch_tris[p.getID()] = []
        for t in range(mesh.ntris):
            patch = self.mesh.getTriPatch(t)
            if patch != None:
                self.patch_tris[patch.getID()].append(t)
                
    def changeColor(self, spec_id, color_idx):
        self.colors[spec] = color_idx
    
    def run(self, end_time, update_interval, record_file = None):
        if record_file != None:
            rc_file = open(record_file, 'w')
            rc_file.seek(os.SEEK_END)
        while self.sim.getTime() < end_time:
            cubit.cmd("Graphics Pause")
            cubit.cmd("delete vertex all")
            if record_file != None:
                rc_file.write("Graphics Pause\n")
                rc_file.write("delete vertex all\n")
            self.sim.advance(update_interval)
            for comp in self.comp_tets.values():
                for t in comp:
                    for spec in self.specs:
                        if(self.sim.getTetCount(t, spec) > 0):
                            center = self.mesh.getTetBarycenter(t)
                            scaled_center = [center[0]/self.scale, 
                                            center[1]/self.scale, 
                                            center[2]/self.scale]
                            cmd = "create vertex %f %f %f color id %i" % (scaled_center[0],
                            scaled_center[1], scaled_center[2], self.colors[spec])
                            cubit.cmd(cmd)
                            if record_file != None:
                                rc_file.write(cmd + "\n")

            for patch in self.patch_tris.values():
                for t in patch:
                    for spec in self.specs:
                        if(self.sim.getTriCount(t, spec) > 0):
                            center = self.mesh.getTriBarycenter(t)
                            scaled_center = [center[0]/self.scale, 
                                            center[1]/self.scale, 
                                            center[2]/self.scale]
                            cmd = "create vertex %f %f %f color id %i" % (scaled_center[0],
                            scaled_center[1], scaled_center[2], self.colors[spec])
                            cubit.cmd(cmd)
                            if record_file != None:
                                rc_file.write(cmd + "\n")
            cubit.cmd("display")
            rc_file.write("display\n")
        if record_file != None:
            rc_file.close()
################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
################################################################################            