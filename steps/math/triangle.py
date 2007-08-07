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
    # If we get a 1-dimensional array, make it 2D.
    if len(t.shape) == 1: t = t.reshape(-1,3)
    
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
        A 2-dimensional array of size N_tri * 3 with barycenters.
    """
    return (p[t[:,0],:] + p[t[:,1],:] + p[t[:,2],:]) / 3.0
    
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def normal(p, t):
    """Compute the normal vector(s) for one or more triangles.
    
    These normal vectors of course have length 1.0. 
    
    PARAMETERS:
        p
            An N_pnt * 3 array of points.
        t
            An N_tri * 3 array of integer indices into p.
    
    RETURNS:
        A 2-dimensional array of size N_tri * 3 with the normal vectors.
    """
    #ntri = t.shape[0]
    #nrms = numpy.empty((ntri,3))
    #nrms[:,0] = ((p[t[:,1],1] - p[t[:,0],1]) * (p[t[:,2],2] - p[t[:,0],2])) \
    #          - ((p[t[:,1],2] - p[t[:,0],2]) * (p[t[:,2],1] - p[t[:,0],1]))
    #nrms[:,1] = ((p[t[:,1],2] - p[t[:,0],2]) * (p[t[:,2],0] - p[t[:,0],0])) \
    #          - ((p[t[:,1],0] - p[t[:,0],0]) * (p[t[:,2],2] - p[t[:,0],2]))
    #nrms[:,2] = ((p[t[:,1],0] - p[t[:,0],0]) * (p[t[:,2],1] - p[t[:,0],1])) \
    #          - ((p[t[:,1],1] - p[t[:,0],1]) * (p[t[:,2],0] - p[t[:,0],0]))
    #norm = numpy.sqrt((nrms * nrms).sum(axis = 1))
    #norm = numpy.c_ [ norm, norm, norm ]
    #return nrms / norm
    tmp1 = p[t[:,1],:] - p[t[:,0],:]
    tmp2 = p[t[:,2],:] - p[t[:,0],:]
    tmp = numpy.cross(tmp1,tmp2)
    norm = numpy.sqrt((tmp * tmp).sum(axis = 1))
    return tmp / numpy.c_ [ norm, norm, norm ]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
