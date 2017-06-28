#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math

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

    return (verts,tets)



#(vx,tx) = make_tri_prism(3,10,1.0)
(vx,tx) = make_tri_prism(1,1,1.0)

print 'double COORDS[][3] = {'
for v in vx:
    print '    {'+', '.join(map(str,v))+'},'
print '};'

print 'unsigned int TETINDICES[][4] = {'
for t in tx:
    print '    {'+', '.join(map(str,t))+'},'
print '};'

