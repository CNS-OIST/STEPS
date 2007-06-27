# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""A variety of auxiliary functions for dealing with tetrahedrons.

Currently, all these methods are implemented in Python and we're relying 
on Numpy for speed. If this would ever turn out to be a bottleneck, we'll 
port them to C/C++. 
"""


import numpy
import numpy.linalg.basic


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def vol(p, t):
    """Compute the volume for one or more tetrahedrons.
    
    Parameters:
        p
            An N_pnt * 3 array of points.
        t
            An N_tet * 4 array of integer indices into p.
    
    Returns:
        A 1-dimensional array of size N_tet with volumes.
    
    Raises:
        ---
    """
    workhorse = numpy.empty((4, 4))
    numtets = t.shape[0]
    vols = numpy.empty((numtets))
    for row in xrange(0, numtets):
        # Row: 0
        ix = t[row,0]
        workhorse[0,0] = p[ix,0]
        workhorse[0,1] = p[ix,1]
        workhorse[0,2] = p[ix,2]
        workhorse[0,3] = 1.0
        # Row: 1
        ix = t[row,1]
        workhorse[1,0] = p[ix,0]
        workhorse[1,1] = p[ix,1]
        workhorse[1,2] = p[ix,2]
        workhorse[1,3] = 1.0
        # Row: 2
        ix = t[row,2]
        workhorse[2,0] = p[ix,0]
        workhorse[2,1] = p[ix,1]
        workhorse[2,2] = p[ix,2]
        workhorse[2,3] = 1.0
        # Row: 3
        ix = t[row,3]
        workhorse[3,0] = p[ix,0]
        workhorse[3,1] = p[ix,1]
        workhorse[3,2] = p[ix,2]
        workhorse[3,3] = 1.0
        # Compute volume.
        vols[ix] = abs(linalg.det(workhorse.transpose()) / 6.0)
    return vols
    

def barycenter(p, t):
    """Compute the barycenter for one or more tetrahedrons.
    """
    return (p[t[:,0],:] + p[t[:,1],:] + p[t[:,2],:] + p[t[:,3],:]) / 4.0


def barycentric(tcp, p):
    """Transform a 3D point into its barycentric coordinates.
    
    Parameters:
        tcp
            Tetrahedron corner points (a 4*3 array).
        p
            The coordinates that should be transformed.
    """
    # Construct and solve linear system satisfying barycentric coordinates.
    a = numpy.empty((3, 3))
    a[0,:] = tcp[1,:] - tcp[0,:]
    a[1,:] = tcp[2,:] - tcp[0,:]
    a[2,:] = tcp[3,:] - tcp[0,:]
    b = numpy.array(p) - tcp[0,:]
    x = numpy.linalg.basic.solve(a, b)

    # The fourth barycentric coordinate is found by subtracting the other
    # three from 1.0. This can lead to serious truncation errors; mediate
    # this by sorting the coordinates from small to large first.
    x.sort()
    
    # Return them.
    return ( x[0], x[1], x[2], 1.0 - x.sum() )


def insideExcl(tcp, p):
    """Test whether a 3D point is inside a tetrahedron (not including the
    tetrahedron's edges, in other words, strictly inside).
    
    Parameters:
        tcp
            Tetrahedron corner points (a 4*3 array).
        p
            The point to test.
    
    Returns:
        A boolean value.
    
    Raises:
        ---
    """
    bc = xformBarycentric(tcp, p)
    if bc[0] <= 0.0: return False
    if bc[1] <= 0.0: return False
    if bc[2] <= 0.0: return False
    if bc[3] <= 0.0: return False
    return True


def insideIncl(tcp, p):
    """Test whether a 3D point is inside a tetrahedron (including the
    tetrahedron's edges).
    
    Parameters:
        tcp
            Tetrahedron corner points (a 4*3 array).
        p
            The point to test.
    
    Returns:
        A boolean value.
    
    Raises:
        ---
    """
    bc = xformBarycentric(tcp, p)
    if bc[0] < 0.0: return False
    if bc[1] < 0.0: return False
    if bc[2] < 0.0: return False
    if bc[3] < 0.0: return False
    return True


def circumradius(p, t):
    """Compute the radius of the circumsphere for one or more tetrahedrons.
    """
    pass


def circumradius2(p, t):
    """Compute the square radius of the circumsphere for one or more 
    tetrahedrons.
    """
    pass


def shortestEdge(p, t):
    """Compute the length of the shortest edge of one or more tetrahedrons.
    """
    pass


def shortestEdge2(p, t):
    """Compute the square length of the shortest edge of one or more
    tetrahedrons.
    """
    pass


def longestEdge(p, t):
    """Compute the length of the longest edge of one or more tetrahedrons.
    """
    pass


def longestEdge2(p, t):
    """Compute the square length of the longest edge of one or more 
    tetrahedrons.
    """
    pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
