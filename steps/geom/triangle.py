# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""A variety of auxiliary functions for dealing with triangles.

Currently, all these methods are implemented in Python and we're relying 
on NumPy for speed. If this would ever turn out to be a bottleneck, we'll 
port them to C/C++. 
"""


import math
import numpy


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def area(p, t):
    """Compute the area for one or more triangles.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tri * 3 array of integer indices into p.
    
    RETURNS:
        A 1-dimensional array of size N_tri with areas.
    """
    p0 = p[t[:,0],:]
    c = numpy.cross(p[t[:,1],:] - p0, p[t[:,2],:] - p0)
    return 0.5 * numpy.array([ math.sqrt(x) for x in (c*c).sum(axis=1) ])
    

def barycenter(p, t):
    """Compute the barycenter for one or more triangles.
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tri * 3 array of integer indices into p.
    
    RETURNS:
        A 1-dimensional array of size N_tri * 3 with barycenters.
    """
    return (p[t[:,0],:] + p[t[:,1],:] + p[t[:,2],:]) / 3.0


def toBarycentric(tcp, p):
    """Transform a 3D point into its barycentric coordinates.
    
    This method (currently) operates on a single point, unlike other 
    methods in this module.
    
    PARAMETERS:
        tcp
            Triangle corner points (a 3*3 array).
        p
            The point that should be transformed.
    
    RETURNS:
        Guess...
    """
    pass
    

def circumradius(p, t):
    """Compute the radius of the circumcircle for one or more triangles.
    """
    pass


def circumradius2(p, t):
    """Compute the square radius of the circumcircle for one or more 
    triangles.
    """
    pass


def shortestEdge(p, t):
    """Compute the length of the shortest edge of one or more triangles.
    """
    pass


def shortestEdge2(p, t):
    """Compute the square length of the shortest edge of one or more
    triangles.
    """
    pass


def longestEdge(p, t):
    """Compute the length of the longest edge of one or more triangles.
    """
    pass


def longestEdge2(p, t):
    """Compute the square length of the longest edge of one or more 
    triangles.
    """
    pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
