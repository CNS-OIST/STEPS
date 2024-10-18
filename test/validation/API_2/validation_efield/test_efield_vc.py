####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

import unittest

import steps.interface

from steps.geom import *
from steps.model import *
from steps.sim import *
from steps.saving import *
from steps.rng import *
from steps.utils import *

import numpy as np
import re

class TestEfieldVC(unittest.TestCase):

    def setUp(self):
        SetVerbosity(0)

    def make_tri_prism(self, k, m, scale):
        # Generate triangular prism with k line segments along each triangular
        # edge and m along vertical extrusion axis. Each line segment is length scale.

        nvert_per_slice = ((k+1)*(k+2))//2
        ntri_per_slice = k*k

        # make vertices
        verts = [(x,y,z) for z in range(0,m+1) for y in range(0,k+1) for x in range(0,k+1-y)]
                    
        # shear and scale ...
        sqrt3o2 = np.sqrt(3.0)/2
        verts = [(x*scale+0.5*y*scale,sqrt3o2*scale*y,scale*z) for (x,y,z) in verts]
       
        # vertex at (x,y,z), x+y <= k.
        def xyz_to_vidx(x,y,z):
            return z*nvert_per_slice + x + y*k + (y*(3-y))/2

        def vidx_to_xyz(v):
            z = v//nvert_per_slice
            x = v%nvert_per_slice
            y = 0
            mx = k+1
            while x>=mx:
               x -= mx
               mx -= 1
               y += 1

            return (x,y,z)

        # ith triangle vertex indices in slice z
        def tri_verts(i,z):
            r = i%k
            s = i//k
            if r+s<k:
                return [xyz_to_vidx(x,y,z) for (x,y) in [(r,s),(r+1,s),(r,s+1)]]
            else:
                r,s = (k-s,k-r)
                return [xyz_to_vidx(x,y,z) for (x,y) in [(r,s),(r-1,s),(r,s-1)]]

        # tets in simple triangular prism between vertices in triangles p and q:
        def tet_verts(p,q,flip):
            vidx = p + q
            vs = [[0,1,2,3],[1,4,3,2],[2,4,5,3]] if not flip else \
                 [[0,1,2,5],[0,4,5,3],[0,1,5,4]]
            return [[vidx[i] for i in tetv] for tetv in vs]

        tets = [tet
                  for z in range(0,m)
                  for t in range(0,ntri_per_slice)
                  for tet in tet_verts(tri_verts(t,z),tri_verts(t,z+1),((t%k)+(t/k))>=k)]

        #print verts
        #print tets

        # flatten vertices, tets lists and make Tetmesh

        mesh = TetMesh.FromData([u for (x,y,z) in verts for u in [x,y,z]],
                             [i for tet in tets for i in tet])

        with mesh:
            # mark faces as ROI
            rois = {name : TriList() for name in ['bottom', 'top', 'face1', 'face2', 'face3']}
            for tri in mesh.tris:
                xyzs = [vidx_to_xyz(v.idx) for v in tri.verts]

                if all([z==0 for (x,y,z) in xyzs]):
                    rois['bottom'].append(tri)
                if all([z==m for (x,y,z) in xyzs]):
                    rois['top'].append(tri)
                if all([x==0 for (x,y,z) in xyzs]):
                    rois['face1'].append(tri)
                if all([y==0 for (x,y,z) in xyzs]):
                    rois['face2'].append(tri)
                if all([x+y==k for (x,y,z) in xyzs]):
                    rois['face3'].append(tri)

            for (roi, roi_tris) in rois.items():
                ROI(roi_tris, name=roi)

            # mark horizontal slices of vertices as ROI too.
            for s in range(0,m+1):
                ROI(VertList(range(s*nvert_per_slice,(s+1)*nvert_per_slice)), name='slice'+str(s))

        return mesh


    def test_efield_vc(self):
        mesh = self.make_tri_prism(3,10,1.0)
        with mesh:
            interior = Compartment.Create(mesh.tets)
            patch = Patch.Create(mesh.surface, interior)
            memb = Membrane.Create([patch], opt_method=1)

        model = Model()
        with model:
            # need at minimum one species
            sA = Species.Create()

        rng = RNG('r123', 512, 123456789)
        sim = Simulation('Tetexact', model, mesh, rng, calcMembPot=True)

        N = 11
        EF_dt = 1e-6 # 1 microsecond 

        # get means across all horizontal slices
        rs = ResultSelector(sim)

        slices = [(int(roi.name[5:]), roi) for roi in mesh.ALL(ROI) if roi.name.startswith('slice')]
        slices.sort(key=lambda x: x[0])

        result = None
        for _, roi in slices:
            res = rs.SUM(rs.VERTS(roi.verts).V) / len(roi.verts)
            result = res if result is None else result << res

        sim.toSave(result, dt=EF_dt)
       
        sim.newRun()
        sim.EfieldDT = EF_dt

        # initialise potential of mesh vertices to 0V
        sim.memb.Potential = 0.0

        # membrane capacitance
        sim.memb.Capac = 0.0

        # volume resistivity
        sim.memb.VolRes = 1.0

        # clamp top and bottom of prism
        vtop = VertList((vi for tri in mesh.top.tris for vi in tri.verts), mesh=mesh)
        vbottom = VertList((vi for tri in mesh.bottom.tris for vi in tri.verts), mesh=mesh)

        sim.VERTS(vtop).V = 10.0
        sim.VERTS(vtop).VClamped = True

        sim.VERTS(vbottom).V = 0.0
        sim.VERTS(vbottom).VClamped = True

        sim.run(N * EF_dt)

        for vertVs in result.data[0, 2:, :]:
            for i, v in enumerate(vertVs):
                assert (abs(v-i) < 1e-10)


def suite():
    all_tests = []
    all_tests.append(unittest.TestLoader().loadTestsFromTestCase(TestEfieldVC))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
