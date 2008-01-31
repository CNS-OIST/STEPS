# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# This file is part of STEPS.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import math
import numpy
import scipy.io.mio as mio
import scipy.sparse as sparse
import steps.thirdparty.qhull as qhull
import sys


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Sphere(object):

    def __init__(self, r):
        self._r = r
    
    def sigdist(self, p):
        npoints = p.shape[0]
        d = numpy.zeros(npoints)
        for i in xrange(npoints):
            d[i] = math.sqrt(p[i,0]*p[i,0] \
                           + p[i,1]*p[i,1] \
                           + p[i,2]*p[i,2]) \
                 - self._r
        return d

    def relsize(self, p):
        npoints = p.shape[0]
        h = numpy.ones(npoints)
        return h


class Cylinder(object):
    
    def __init__(self, p1, p2, rad):
        self._p1 = numpy.asarray(p1)
        self._p2 = numpy.asarray(p2)
        self._rad = rad
        self._p12 = self._p2 - self._p1
        self._dp2 = numpy.sum(self._p12 * self._p12)
        self._dp = math.sqrt(self._dp2)
        
    def sigdist(self, p):
        """Returns the signed distance of one or more test points.
        """
        npoints = p.shape[0] 
        # Compute projected parametric coordinates.
        t = -numpy.sum((self._p1 - p) * self._p12, axis=1) / self._dp2
        # Compute coaxial signed distance to p1 and p2.
        d0 = -t * self._dp
        d1 = (t - 1.0) * self._dp
        # Compute radial distance.
        p00 = (numpy.c_[t,t,t] * self._p12) + self._p1
        pp00 = p00-p
        dr = numpy.sqrt(numpy.sum(pp00*pp00,axis=1)) - self._rad
        df = numpy.c_[d0,d1,dr]
        df.sort(axis=1)
        return df[:,2]
    
    def relsize(self, p):
        npoints = p.shape[0]
        h = numpy.ones(npoints)
        return h

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def unique(s):
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    
    Source: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52560
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u


def meshgen(gd, p0, h):
    ptol = 0.001
    ttol = 0.1
    L0mult = 1.1
    deltat = 0.1
    geps = 0.1 * h
    deps = math.sqrt(numpy.finfo(float).eps) * h
    
    # 2. Remove points outside the region (no rejection method).
    p = p0.compress(gd.sigdist(p0) < geps, 0)
    N = p.shape[0]
    
    count = 0
    p0 = numpy.inf
    while True:
        # 3. Retriangulation by Delaunay.
        if numpy.max(numpy.sqrt(numpy.sum((p-p0)**2,1))) > ttol * h:
            p0 = p.copy()
            # NOTE: get rid of the delny stuff, we only need 't'. Everything
            # else is a waste of time. 
            # UPDATE: we actually also need a list of neighbours, as in step 4.
            (n,f,t) = qhull.delny(p)
            print '====='
            t = numpy.array(t)
            print t.shape
            t.astype(int)
            pmid = numpy.zeros((t.shape[0],3))
            pmid = pmid + p[t[:,0],:] + p[t[:,1],:] + p[t[:,2],:] + p[t[:,3],:]
            pmid = pmid / 4.0
            # Remove outside tetrahedrons.
            print t.shape
            t = t.compress(gd.sigdist(pmid) < -geps, 0)
            print t.shape
            
            # 4. Describe each edge by a unique pair of nodes.
            # CLEAN UP!
            pair = numpy.zeros((0,2),dtype=int)
            localpairs = numpy.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
            for ii in xrange(localpairs.shape[0]):
                pair = numpy.r_[pair, t[:,localpairs[ii,:]]]
            pair2 = numpy.sort(pair, 1).tolist()
            pair2 = unique(pair2)
            pair = numpy.array(pair2)
            
            # 5. Graphical output of the current mesh.
            # SKIP
            
            # Increase triangulation counter.
            count = count + 1
        
        print t.shape
        
        # 6. Move mesh points based on edge lengths L and forces F.
        #
        # TODO: this whole section would probably benefit from a clean
        # up, as it's currently a 'literal' translation of the original
        # Matlab code. Some important Matlab functions are not available
        # in NumPy/SciPy however, which explain the awkwardness of the code
        # below. See the annotated printout of 'distmeshnd.m' for notes.
        bars = p[pair[:,0],:] - p[pair[:,1],:]
        L = numpy.sqrt(numpy.sum(bars**2,1))
        L0 = gd.relsize((p[pair[:,0],:] + p[pair[:,1],:]) / 2)
        L0 = L0 * L0mult * (numpy.sum(L**3)/numpy.sum(L0**3)) ** (1.0/3.0)
        F = numpy.clip(L0-L, 0.0, numpy.inf)
        # Fbar = numpy.c_[bars,-bars] * numpy.repmat(F/L, 1, 6)
        # --> NumPy doesn't seem to implement repmat (yet?)
        FL = F/L
        Fbar = numpy.c_[bars,-bars] * numpy.c_[FL, FL, FL, FL, FL, FL]
        tmp_i = pair[:,[0,0,0,1,1,1]]
        # tmp_j = numpy.repmat(numpy.array([0,1,2,0,1,2]), pair.shape[0], 1)
        # --> NumPy doesn't seem to implement repmat (yet?)
        tmp_j = numpy.repeat(numpy.array([[0,1,2,0,1,2]]),pair.shape[0],0)
        ij = numpy.c_[tmp_i.flatten(),tmp_j.flatten()]
        ij = ij.transpose()
        dp = sparse.csc_matrix((Fbar.flatten(),ij), dims=(N,3)).toarray()
        # dp[1,size(fix,1),:]=0
        p = p + deltat * dp
        
        # 7. Bring outside points back to the boundary.
        d = gd.sigdist(p)
        ix = (d > 0.0)
        nix = numpy.sum(ix)
        gradd = numpy.zeros((nix, 3))
        for ii in [0, 1, 2]:
            a = numpy.zeros((nix,3))
            a[:,ii] = deps
            d1x = gd.sigdist(p[ix,:] + a)
            gradd[:,ii] = (d1x-d[ix]) / deps;
        p[ix,:] = p[ix,:] - (numpy.c_[d[ix],d[ix],d[ix]] * gradd)
        
        # 8. Termination criterion.
        maxdp = numpy.max(deltat * numpy.sqrt(numpy.sum(dp[d<-geps]**2,1)))
        print '%f < %f ?' % (maxdp, ptol * h)
        if maxdp < (ptol * h):
            print t.shape
            return (p,t)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def extractSurfTri(p,t):
    """Extract the surface triangles from an N*4 array of indices and orient
    them properly.
    """
    # Check for the trivial case.
    if t.shape[0] == 0:
        return
    
    # Combine all possible faces.
    faces = numpy.r_[ t[:,[0,1,2]], t[:,[0,1,3]], t[:,[0,2,3]], t[:,[1,2,3]] ]
    inods = numpy.r_[ t[:,3], t[:,2], t[:,1], t[:,0] ]
    faces.sort(1)
    
    # 
    tricnt = dict()
    for i in xrange(0, faces.shape[0]):
        curtri = ( faces[i,0], faces[i,1], faces[i,2] )
        n = tricnt.setdefault(curtri, 0)
        tricnt[curtri] = n + 1
    surftriidx = [ ]
    for i in xrange(0, faces.shape[0]):
        curtri = ( faces[i,0], faces[i,1], faces[i,2] )
        if tricnt[curtri] == 1:
            surftriidx.append(i)
    
    # Fetch surface triangles and opposing (interal) nodes.
    faces = faces[surftriidx,:]
    inods = inods[surftriidx]
    
    # Handle the orientation
    v1 = p[faces[:,1],:] - p[faces[:,0],:]
    v2 = p[faces[:,2],:] - p[faces[:,0],:]
    v3 = p[inods,:] - p[faces[:,0],:]
    cr = numpy.cross(v1,v2)
    print cr.shape
    print v3.shape
    dt = numpy.sum(cr * v3,1)
    print dt.shape
    ix = dt > 0
    tmp = faces[ix,2]
    faces[ix,2] = faces[ix,1]
    faces[ix,1] = tmp
    return faces
    

def meshgenloop(maxcount, gd, p0, h):
    ptol = 0.001
    ttol = 0.1
    L0mult = 1.1
    deltat = 0.1
    geps = 0.1 * h
    deps = math.sqrt(numpy.finfo(float).eps) * h
    
    # 2. Remove points outside the region (no rejection method).
    p = p0.compress(gd.sigdist(p0) < geps, 0)
    N = p.shape[0]
    
    count = 0
    p0 = numpy.inf
    while True:
        # 3. Retriangulation by Delaunay.
        if numpy.max(numpy.sqrt(numpy.sum((p-p0)**2,1))) > ttol * h:
            p0 = p.copy()
            # NOTE: get rid of the delny stuff, we only need 't'. Everything
            # else is a waste of time. 
            # UPDATE: we actually also need a list of neighbours, as in step 4.
            (n,f,t) = qhull.delny(p)
            print '====='
            t = numpy.array(t)
            print t.shape
            t.astype(int)
            pmid = numpy.zeros((t.shape[0],3))
            pmid = pmid + p[t[:,0],:] + p[t[:,1],:] + p[t[:,2],:] + p[t[:,3],:]
            pmid = pmid / 4.0
            # Remove outside tetrahedrons.
            print t.shape
            t = t.compress(gd.sigdist(pmid) < -geps, 0)
            print t.shape
            
            # 4. Describe each edge by a unique pair of nodes.
            # CLEAN UP!
            pair = numpy.zeros((0,2),dtype=int)
            localpairs = numpy.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
            for ii in xrange(localpairs.shape[0]):
                pair = numpy.r_[pair, t[:,localpairs[ii,:]]]
            pair2 = numpy.sort(pair, 1).tolist()
            pair2 = unique(pair2)
            pair = numpy.array(pair2)
            
            # 5. Graphical output of the current mesh.
            # SKIP
            
            # Increase triangulation counter.
            if count == maxcount:
                return (p, t)
            count = count + 1
        
        print t.shape
        
        # 6. Move mesh points based on edge lengths L and forces F.
        #
        # TODO: this whole section would probably benefit from a clean
        # up, as it's currently a 'literal' translation of the original
        # Matlab code. Some important Matlab functions are not available
        # in NumPy/SciPy however, which explain the awkwardness of the code
        # below. See the annotated printout of 'distmeshnd.m' for notes.
        bars = p[pair[:,0],:] - p[pair[:,1],:]
        L = numpy.sqrt(numpy.sum(bars**2,1))
        L0 = gd.relsize((p[pair[:,0],:] + p[pair[:,1],:]) / 2)
        L0 = L0 * L0mult * (numpy.sum(L**3)/numpy.sum(L0**3)) ** (1.0/3.0)
        F = numpy.clip(L0-L, 0.0, numpy.inf)
        # Fbar = numpy.c_[bars,-bars] * numpy.repmat(F/L, 1, 6)
        # --> NumPy doesn't seem to implement repmat (yet?)
        FL = F/L
        Fbar = numpy.c_[bars,-bars] * numpy.c_[FL, FL, FL, FL, FL, FL]
        tmp_i = pair[:,[0,0,0,1,1,1]]
        # tmp_j = numpy.repmat(numpy.array([0,1,2,0,1,2]), pair.shape[0], 1)
        # --> NumPy doesn't seem to implement repmat (yet?)
        tmp_j = numpy.repeat(numpy.array([[0,1,2,0,1,2]]),pair.shape[0],0)
        ij = numpy.c_[tmp_i.flatten(),tmp_j.flatten()]
        ij = ij.transpose()
        dp = sparse.csc_matrix((Fbar.flatten(),ij), dims=(N,3)).toarray()
        # dp[1,size(fix,1),:]=0
        p = p + deltat * dp
        
        # 7. Bring outside points back to the boundary.
        d = gd.sigdist(p)
        ix = (d > 0.0)
        nix = numpy.sum(ix)
        gradd = numpy.zeros((nix, 3))
        for ii in [0, 1, 2]:
            a = numpy.zeros((nix,3))
            a[:,ii] = deps
            d1x = gd.sigdist(p[ix,:] + a)
            gradd[:,ii] = (d1x-d[ix]) / deps;
        p[ix,:] = p[ix,:] - (numpy.c_[d[ix],d[ix],d[ix]] * gradd)
        
        # 8. Termination criterion.
        maxdp = numpy.max(deltat * numpy.sqrt(numpy.sum(dp[d<-geps]**2,1)))
        print '%f < %f ?' % (maxdp, ptol * h)
        if maxdp < (ptol * h):
            print t.shape
            return (p,t)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


_vis_pnts = []
_vis_tets = []
_vis_tets2 = []
_vis_tris = []

_vis_lastx = 0
_vis_lasty = 0

_vis_ccnt = 0


def norm3(p1,p2,p3):
    v1 = (p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
    v2 = (p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2])
    return ((v1[1]*v2[2]) - (v1[2]*v2[1]), \
            (v1[2]*v2[0]) - (v1[0]*v2[2]), \
            (v1[0]*v2[1]) - (v1[1]*v2[0]))


def handleGLUTDisplay():
    global _vis_pnts, _vis_tets, _vis_tets2, _vis_tris
    global _vis_lastx, _vis_lasty
    
    # Clear screen.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glColor3f(1.0, 1.0, 1.0)
    
    # Set up drawing.
    glLoadIdentity()
    gluLookAt(0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    glRotatef(_vis_lastx, 0.0, 1.0, 0.0)
    glRotatef(_vis_lasty, 1.0, 0.0, 0.0)
    
    # Plot vertex points.
    glPointSize(2.0)
    glBegin(GL_POINTS)
    for i in xrange(0, _vis_pnts.shape[0]):
        glVertex3f(_vis_pnts[i,0], _vis_pnts[i,1], _vis_pnts[i,2])
    glEnd()
    
    # Plot tetrahedrons points.
    glPointSize(1.0)
    glColor3f(1.0, 1.0, 1.0)
    glBegin(GL_TRIANGLES)
    for i in xrange(0, _vis_tris.shape[0]):
        t0 = _vis_tris[i,0]
        t1 = _vis_tris[i,1]
        t2 = _vis_tris[i,2]
        p0 = (_vis_pnts[t0,0],_vis_pnts[t0,1],_vis_pnts[t0,2])
        p1 = (_vis_pnts[t1,0],_vis_pnts[t1,1],_vis_pnts[t1,2])
        p2 = (_vis_pnts[t2,0],_vis_pnts[t2,1],_vis_pnts[t2,2])
        glNormal3fv(norm3(p0,p1,p2))
        glVertex3fv(p0)
        glVertex3fv(p1)
        glVertex3fv(p2)
    glEnd()
    
    # Draw.
    glutSwapBuffers()


def handleGLUTReshape(w, h):
    global _vis_lastx, _vis_lasty
    glViewport(0, 0, w, h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0)
    glMatrixMode(GL_MODELVIEW)


def handleGLUTKeyboard(key, x, y):
    global _vis_pnts, _vis_tets, _vis_tets2, _vis_tris
    global _vis_ccnt
    
    if key=='s':
        (bli,bla,blo) = qhull.delny(_vis_pnts)
        _vis_tets = numpy.array(blo)
        dic = { 't':_vis_tets, 'p':_vis_pnts }
        mio.savemat('/Users/stefan/mesh.mat', dic)
        return
    
    # Refine the mesh.
    g = Sphere(1.0)
    p0 = _vis_pnts.copy()
    ( _vis_pnts, _vis_tets ) = meshgenloop(1, g, p0, 0.2)
    (bli,bla,blo) = qhull.delny(_vis_pnts)
    
    # What'll be rendered?
    _vis_tets = numpy.array(blo)
    _vis_tets2 = []
    for i in xrange(0, _vis_tets.shape[0]):
        t0 = _vis_tets[i,0]
        t1 = _vis_tets[i,1]
        t2 = _vis_tets[i,2]
        t3 = _vis_tets[i,3]
        p0 = (_vis_pnts[t0,0],_vis_pnts[t0,1],_vis_pnts[t0,2])
        if p0[2] > 0.4: continue
        p1 = (_vis_pnts[t1,0],_vis_pnts[t1,1],_vis_pnts[t1,2])
        if p1[2] > 0.4: continue
        p2 = (_vis_pnts[t2,0],_vis_pnts[t2,1],_vis_pnts[t2,2])
        if p2[2] > 0.4: continue
        p3 = (_vis_pnts[t3,0],_vis_pnts[t3,1],_vis_pnts[t3,2])
        if p3[2] > 0.4: continue
        _vis_tets2.append([t0,t1,t2,t3])
    _vis_tets2 = numpy.array(_vis_tets2)
    _vis_tris = extractSurfTri(_vis_pnts,_vis_tets2)
    _vis_tris = numpy.array(_vis_tris)
    
    # Set window title.
    _vis_ccnt = _vis_ccnt + 1
    glutSetWindowTitle("Resulting mesh (iteration %d)" % (_vis_ccnt))
    
    # Request a redisplay.
    glutPostRedisplay()


def handleGLUTMouseMotion(x, y):
    global _vis_lastx, _vis_lasty
    _vis_lastx = x
    _vis_lasty = y
    glutPostRedisplay()


def runtest(n):
    global _vis_pnts, _vis_tets, _vis_tets2, _vis_tris
    global _vis_ccnt
    
    # Generate the mesh.
    g = Sphere(1.0)
    _vis_pnts = (numpy.random.rand(n,3) * 2.0) - 1.0
    ( _vis_pnts, _vis_tets ) = meshgenloop(1, g, _vis_pnts, 0.2)
    _vis_tets2 = _vis_tets.copy()
    _vis_tris = extractSurfTri(_vis_pnts,_vis_tets2)
    _vis_tris = numpy.array(_vis_tris)
    _vis_ccnt = 1
    
    # Intialize GLUT
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(800, 800)
    glutCreateWindow("Resulting mesh (iteration %d)" % (_vis_ccnt))
    glutDisplayFunc(handleGLUTDisplay)
    glutKeyboardFunc(handleGLUTKeyboard)
    glutReshapeFunc(handleGLUTReshape)
    glutMotionFunc(handleGLUTMouseMotion)
    
    # Initialize GL
    #glClearDepth(1.0)
    glEnable(GL_DEPTH_TEST)
    glClearColor(0.0, 0.0, 0.0, 0.0)
    glShadeModel(GL_FLAT)
    
    # Set view parameters.
    #glMatrixMode(GL_PROJECTION)
    #glLoadIdentity()
    #glFrustum(-9.0, 9.0, -9.0, 9.0, 50.0, 150.0)
    #glMatrixMode(GL_MODELVIEW)
    
    # Set Lighting & elementary shading.
    glLightfv(GL_LIGHT0, GL_POSITION, (0.0, 0.0, 50.0, 0.0))
    glLightfv(GL_LIGHT0, GL_AMBIENT, (0.4, 0.4, 0.4, 1.0))
    glLightfv(GL_LIGHT0, GL_DIFFUSE, (1.0, 1.0, 1.0, 1.0))
    glEnable(GL_LIGHT0)
    glEnable(GL_LIGHTING)
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE)
    glEnable(GL_COLOR_MATERIAL)
    glEnable(GL_CULL_FACE)
    glCullFace(GL_BACK)
    
    # Run the visualization loop.
    glutMainLoop()
    
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
