# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""A variety of auxiliary functions for dealing with 3D triangles.

The following sets of functionality are offered by this module:
    - Computing triangle areas
    - Barycentric centroids

Currently, all these methods are implemented in Python and we're relying 
on NumPy for speed. If this would ever turn out to be a bottleneck, we'll 
port (parts of) them to C/C++.

SEE ALSO:
    steps.math.tetrahedron
    steps.math.triangle_test
"""


import numpy
import numpy.linalg as linalg


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
    return 0.5 * numpy.sqrt((c * c).sum(axis = 1))
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
