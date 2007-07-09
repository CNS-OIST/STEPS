# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""A variety of auxiliary functions for dealing with tetrahedrons.

Currently, all these methods are implemented in Python and we're relying 
on NumPy for speed. If this would ever turn out to be a bottleneck, we'll 
port them to C/C++. 
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
