#!/usr/bin/env python
# -*- coding: utf-8 -*-

import steps.quiet
import steps.geom as sgeom
import steps.model as smodel
import steps.solver as ssolver
import steps.rng as srng
import numpy as np
import math
import re

def make_tri_prism(k,m,scale):
    # Generate triangular prism with k line segments along each triangular
    # edge and m along vertical extrusion axis. Each line segment is length scale.

    nvert_per_slice = ((k+1)*(k+2))/2
    ntri_per_slice = k*k

    # make vertices
    verts = [(x,y,z) for z in range(0,m+1) for y in range(0,k+1) for x in range(0,k+1-y)]
                
    # shear and scale ...
    sqrt3o2 = math.sqrt(3.0)/2
    verts = [(x*scale+0.5*y*scale,sqrt3o2*scale*y,scale*z) for (x,y,z) in verts]
   
    # vertex at (x,y,z), x+y <= k.
    def xyz_to_vidx(x,y,z):
        return z*nvert_per_slice + x + y*k + (y*(3-y))/2

    def vidx_to_xyz(v):
        z = v/nvert_per_slice
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
        s = i/k
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

    mesh = sgeom.Tetmesh([u for (x,y,z) in verts for u in [x,y,z]],
                         [i for tet in tets for i in tet])

    # mark faces as ROI
    rois = { 'bottom': [], 'top': [], 'face1': [], 'face2': [], 'face3': []}
    for t in xrange(0,mesh.ntris):
        xyzs = [vidx_to_xyz(v) for v in mesh.getTri(t)]

        if all([z==0 for (x,y,z) in xyzs]): rois['bottom'].append(t)
        if all([z==m for (x,y,z) in xyzs]): rois['top'].append(t)
        if all([x==0 for (x,y,z) in xyzs]): rois['face1'].append(t)
        if all([y==0 for (x,y,z) in xyzs]): rois['face2'].append(t)
        if all([x+y==k for (x,y,z) in xyzs]): rois['face3'].append(t)

    for (roi,roi_tris) in rois.iteritems():
        mesh.addROI(roi,sgeom.ELEM_TRI,roi_tris)

    # mark horizontal slices of vertices as ROI too.
    for s in range(0,m+1):
        mesh.addROI('slice'+str(s),sgeom.ELEM_VERTEX,range(s*nvert_per_slice,(s+1)*nvert_per_slice))

    return mesh


def test_efield_vc():
    mesh = make_tri_prism(3,10,1.0)
    interior = sgeom.TmComp('interior',mesh,range(mesh.ntets))

    model = smodel.Model()

    # need at minimum one species
    sA = smodel.Spec('A',model)

    patch = sgeom.TmPatch('patch',mesh,mesh.getSurfTris(),interior)
    memb = sgeom.Memb('membrane',mesh,[patch],opt_method=1)

    rng = srng.create('r123',512)
    sim = ssolver.Tetexact(model,mesh,rng,True)

    EF_dt = 1e-6 # 1 microsecond 
   
    sim.reset()
    sim.setEfieldDT(EF_dt)

    # initialise potential of mesh vertices to 0V
    sim.setMembPotential('membrane',0)

    # membrane capacitance
    sim.setMembCapac('membrane',0)

    # volume resistivity
    sim.setMembVolRes('membrane',1) # 1 ohm·m

    # clamp top and bottom of prism
    vtop = set((vi for tri in mesh.getROIData('top') for vi in mesh.getTri(tri)))
    vbottom = set((vi for tri in mesh.getROIData('bottom') for vi in mesh.getTri(tri)))

    for v in vtop:
        sim.setVertV(v,10.0)
        sim.setVertVClamped(v,True)

    for v in vbottom:
        sim.setVertV(v,0.0)
        sim.setVertVClamped(v,True)

    # get means across all horizontal slices
    slices = dict((int(mo.group(1)),mesh.getROIData(s)) for s in mesh.getAllROINames() for mo in [re.search('slice(\d+)',s)] if mo)
    slice_keys = slices.keys()
    slice_keys.sort()

    N = 10 
    result = np.zeros((N,1+len(slice_keys)))

    for k in xrange(N):
       result[k,0] = k
       sim.advance(EF_dt);
       j = 1
       for s in slice_keys:
           nv = 0
           sumv = 0.0
           for vi in slices[s]:
               nv += 1
               sumv += sim.getVertV(vi)

           result[k,j]=sumv/nv
           j += 1

    for i in np.arange(1.0, 12.0):
        for j in result[1:,int(i)] : assert (abs(j-(i-1)) < 1e-10 )


test_efield_vc()

