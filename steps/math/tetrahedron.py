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


"""A variety of auxiliary functions for dealing with tetrahedrons.

The following sets of functionality are offered by this module:
    - Computing tetrahedron volumes
    - Barycentric coordinates
    - Circumspheres
    - Edge lengths

Currently, all these methods are implemented in Python and we're relying 
on NumPy for speed. If this would ever turn out to be a bottleneck, we'll 
port (parts of) them to C/C++.

SEE ALSO:
    steps.math.tetrahedron_test
    steps.math.triangle
"""


import numpy
import numpy.linalg as linalg


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def vol(p, t):
    """Compute the volume for one or more tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1-dimensional array of size N_tet with volumes.
    """
    # If we get a 1-dimensional array, make it 2D.
    if len(t.shape) == 1: t = t.reshape(-1,4)
    
    workhorse = numpy.empty((4, 4))
    numtets = t.shape[0]
    vols = numpy.empty((numtets))
    for ix in xrange(0, numtets):
        # Fill in the rows with points.
        workhorse[:,0:3] = p[t[ix,:],:]
        # Set final column to 1.0
        workhorse[:,3] = 1.0
        # Compute volume.
        vols[ix] = abs(linalg.det(workhorse.transpose()) / 6.0)
    return vols
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def barycenter(p, t):
    """Compute the barycenter for one or more tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1-dimensional array of size N_tet * 3 with barycenters.
    """
    # If we get a 1-dimensional array, make it 2D.
    if len(t.shape) == 1: t = t.reshape(-1,4)
    
    return (p[t[:,0],:] + p[t[:,1],:] + p[t[:,2],:] + p[t[:,3],:]) / 4.0


def toBarycentric(tcp, p):
    """Transform a 3D point into its barycentric coordinates, determined
    by the corner points of a tetrahedron.
    
    Unlike other methods in this module, this method (currently) operates 
    on a single point.
    
    PARAMETERS:
        tcp
            Tetrahedron corner points (a 4*3 array).
        p
            The point that should be transformed.
    
    RETURNS:
        Bary centric coordinates.
    
    TODO:
        The precision of this method appears to be lacking a bit at this 
        moment, at least when comparing to results obtained with Matlab.
        (In tetrahedron_test.py, we can only expect tolerance to 3 digits
        when comparing the ratio of STEPS- and Matlab computed values
        to 1.0...!) The sort() doesn't seem to be doing much about this, 
        so we'll have to look into what linalg.solve() does, and how to 
        improve on that.
    """
    # Construct and solve linear system satisfying barycentric coordinates.
    a = numpy.empty((3, 3))
    a[:,0] = (tcp[1,:] - tcp[0,:]).transpose()
    a[:,1] = (tcp[2,:] - tcp[0,:]).transpose()
    a[:,2] = (tcp[3,:] - tcp[0,:]).transpose()
    b = numpy.array(p) - tcp[0,:]
    x = linalg.solve(a, b)

    # The fourth barycentric coordinate is found by subtracting the other
    # three from 1.0. This can lead to serious truncation errors; mediate
    # this by sorting the coordinates from small to large first.
    #x.sort()
    
    # Return them.
    return ( x[0], x[1], x[2], 1.0 - x.sum() )


def inside(tcp, p):
    """Test whether a 3D point is inside a tetrahedron.
    
    This method (currently) operates on a single point, unlike other 
    methods in this module.
    
    PARAMETERS:
        tcp
            Tetrahedron corner points (a 4*3 array).
        p
            The point to test.
    
    RETURNS:
        A boolean value.
    
    RAISES:
        ---
    """
    bc = toBarycentric(tcp, p)
    if bc[0] < 0.0: return False
    if bc[1] < 0.0: return False
    if bc[2] < 0.0: return False
    if bc[3] < 0.0: return False
    return True


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def circumsphere(p, t):
    """Compute the circumsphere for one or more tetrahedrons.
    
    The circumsphere of a tetrahedron passes through all 4 corner points.
    It can be found by solving the linear system A.X = B, with A given by:
    
        |  x2 - x1   x3 - x1   x4 - x1  |
        |  y2 - y1   y3 - y1   y4 - y1  |
        |  z2 - z1   z3 - z1   z4 - z1  |
    
    and B:
    
        | (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 |
        | (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 |
        | (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 |
    
    The radius is computed using the solution (x,y,z) to the linear system:
    
        r = 1/2 * sqrt(x^2 + y^2 + z^2)
        
    The coordinates of the sphere's center:
    
            | x1 |         | x |
        C = | y1 | + 1/2 * | y |
            | z1 |         | z |
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A N_tet * 4 array of floating point values. Columns 0 to 2
        are the coordinates of the circumsphere's centers, column 3
        are the radii.
    
    NOTES:
        Even though this is the standard way of computing the properties
        of a circumsphere, imo this has one main problem: it involves a 
        lot of subtraction, at least in the part where the linear system 
        is being set up. This might have an impact on the precision of the 
        result. Maybe there are better ways of doing it? (Since this method
        will be mostly used in quality measures, i.e. outside of inner 
        loops or responsive GUI code, efficiency is of less importance.)
        
        In general: results appear to become worse with increased 
        coplanarity.
    
    SEE ALSO:
        steps.math.tetrahedron.circumradius()
        steps.math.tetrahedron.circumradius2()
    """
    ntets = t.shape[0]
    sol = numpy.empty((ntets,3))
    for ii in xrange(0, ntets):
        a = p[t[ii,1:4],:]
        t1 = p[t[ii,0],:]
        a[0,:] = a[0,:] - t1
        a[1,:] = a[1,:] - t1
        a[2,:] = a[2,:] - t1
        b = (a * a).sum(axis = 1)
        sol[ii,:] = linalg.solve(a, b)
    res = numpy.empty((ntets,4))
    res[:,0:3] = p[t[:,0],:] + (0.5 * sol)
    res[:,3] = 0.5 * numpy.sqrt((sol * sol).sum(axis = 1))
    return res


def circumradius(p, t):
    """Compute the radius of the circumsphere for one or more tetrahedrons.
    
    This value can be used for computing tetmesh quality measures.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1D array of N_tet radii.
    
    SEE ALSO:
        steps.math.tetrahedron.circumsphere()
        steps.math.tetrahedron.circumradius2()
    """
    ntets = t.shape[0]
    sol = numpy.empty((ntets,3))
    for ii in xrange(0, ntets):
        a = p[t[ii,1:4],:]
        t1 = p[t[ii,0],:]
        a[0,:] = a[0,:] - t1
        a[1,:] = a[1,:] - t1
        a[2,:] = a[2,:] - t1
        b = (a * a).sum(axis = 1)
        sol[ii,:] = linalg.solve(a, b)
    return 0.5 * numpy.sqrt((sol * sol).sum(axis = 1))


def circumradius2(p, t):
    """Compute the square radius of the circumsphere for one or more 
    tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1D array of N_tet square radii.
    
    SEE ALSO:
        steps.math.tetrahedron.circumsphere()
        steps.math.tetrahedron.circumradius()
    """
    ntets = t.shape[0]
    sol = numpy.empty((ntets,3))
    for ii in xrange(0, ntets):
        a = p[t[ii,1:4],:]
        t1 = p[t[ii,0],:]
        a[0,:] = a[0,:] - t1
        a[1,:] = a[1,:] - t1
        a[2,:] = a[2,:] - t1
        b = (a * a).sum(axis = 1)
        sol[ii,:] = linalg.solve(a, b)
    return 0.25 * (sol * sol).sum(axis = 1)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def shortestEdge(p, t):
    """Compute the length of the shortest edge of one or more tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1D array of N_tet shortest edge lengths.
    
    SEE ALSO:
        steps.math.tetrahedron.shortestEdge2()
        steps.math.tetrahedron.longestEdge()
        steps.math.tetrahedron.longestEdge2()
    """
    ntets = t.shape[0]
    edges = numpy.empty((6,ntets))
    edges[0,:] = numpy.square(p[t[:,1],:] - p[t[:,0],:]).sum(axis = 1)
    edges[1,:] = numpy.square(p[t[:,2],:] - p[t[:,0],:]).sum(axis = 1)
    edges[2,:] = numpy.square(p[t[:,2],:] - p[t[:,1],:]).sum(axis = 1)
    edges[3,:] = numpy.square(p[t[:,3],:] - p[t[:,0],:]).sum(axis = 1)
    edges[4,:] = numpy.square(p[t[:,3],:] - p[t[:,1],:]).sum(axis = 1)
    edges[5,:] = numpy.square(p[t[:,3],:] - p[t[:,2],:]).sum(axis = 1)
    return numpy.sqrt(edges.min(axis=0))


def shortestEdge2(p, t):
    """Compute the square length of the shortest edge of one or more
    tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1D array of N_tet shortest edge square lengths.
    
    SEE ALSO:
        steps.math.tetrahedron.shortestEdge()
        steps.math.tetrahedron.longestEdge()
        steps.math.tetrahedron.longestEdge2()
    """
    ntets = t.shape[0]
    edges = numpy.empty((6,ntets))
    edges[0,:] = numpy.square(p[t[:,1],:] - p[t[:,0],:]).sum(axis = 1)
    edges[1,:] = numpy.square(p[t[:,2],:] - p[t[:,0],:]).sum(axis = 1)
    edges[2,:] = numpy.square(p[t[:,2],:] - p[t[:,1],:]).sum(axis = 1)
    edges[3,:] = numpy.square(p[t[:,3],:] - p[t[:,0],:]).sum(axis = 1)
    edges[4,:] = numpy.square(p[t[:,3],:] - p[t[:,1],:]).sum(axis = 1)
    edges[5,:] = numpy.square(p[t[:,3],:] - p[t[:,2],:]).sum(axis = 1)
    return edges.min(axis=0)


def longestEdge(p, t):
    """Compute the length of the longest edge of one or more tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1D array of N_tet longest edge lengths.
    
    SEE ALSO:
        steps.math.tetrahedron.shortestEdge()
        steps.math.tetrahedron.shortestEdge2()
        steps.math.tetrahedron.longestEdge2()
    """
    ntets = t.shape[0]
    edges = numpy.empty((6,ntets))
    edges[0,:] = numpy.square(p[t[:,1],:] - p[t[:,0],:]).sum(axis = 1)
    edges[1,:] = numpy.square(p[t[:,2],:] - p[t[:,0],:]).sum(axis = 1)
    edges[2,:] = numpy.square(p[t[:,2],:] - p[t[:,1],:]).sum(axis = 1)
    edges[3,:] = numpy.square(p[t[:,3],:] - p[t[:,0],:]).sum(axis = 1)
    edges[4,:] = numpy.square(p[t[:,3],:] - p[t[:,1],:]).sum(axis = 1)
    edges[5,:] = numpy.square(p[t[:,3],:] - p[t[:,2],:]).sum(axis = 1)
    return numpy.sqrt(edges.max(axis=0))


def longestEdge2(p, t):
    """Compute the square length of the longest edge of one or more 
    tetrahedrons.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    RETURNS:
        A 1D array of N_tet longest square edge lengths.
    
    SEE ALSO:
        steps.math.tetrahedron.shortestEdge()
        steps.math.tetrahedron.shortestEdge2()
        steps.math.tetrahedron.longestEdge()
    """
    ntets = t.shape[0]
    edges = numpy.empty((6,ntets))
    edges[0,:] = numpy.square(p[t[:,1],:] - p[t[:,0],:]).sum(axis = 1)
    edges[1,:] = numpy.square(p[t[:,2],:] - p[t[:,0],:]).sum(axis = 1)
    edges[2,:] = numpy.square(p[t[:,2],:] - p[t[:,1],:]).sum(axis = 1)
    edges[3,:] = numpy.square(p[t[:,3],:] - p[t[:,0],:]).sum(axis = 1)
    edges[4,:] = numpy.square(p[t[:,3],:] - p[t[:,1],:]).sum(axis = 1)
    edges[5,:] = numpy.square(p[t[:,3],:] - p[t[:,2],:]).sum(axis = 1)
    return edges.max(axis=0)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
