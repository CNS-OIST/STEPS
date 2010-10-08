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
    tet_list = []
    for t in range(ntets):
        center = mesh.getTetBarycenter(t)
        cubit_center = [center[0]/scale, center[1]/scale, center[2]/scale]
        status = body.point_containment(cubit_center)
        if status == 1 or status == 2:
            tet_list.append(t)
    return tet_list
        
def highlightTets(steps_tets, tet_proxy):
    cmd_str = "highlight tet "
    for st in steps_tets:
        cubit_id = tet_proxy.getImportID(st)
        cubit.cmd("highlight tet %i" % (cubit_id))