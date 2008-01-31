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


from pylab import *
from steps.error import ArgumentError

import math
import numpy


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def makeContour(p0, uaxis, vaxis, uvals, vvals, func):
    """Sample a 3D function over a cross section plane. This plane is 
    defined by its origin (p0) and two 3D increment vectors dx and dy.
    The resulting data can be used in Matplotlib's contour() function
    (hence the name).
    
    PARAMETERS:
        p0
            The origin of the planar cross section. Must be a 3D point,
            either in the form of a set, list or Numpy array.
        uaxis, vaxis
            The primary axis vectors of the plane. Both are 3D vectors,
            specified either as a set, list or NumPy array. They don't
            have to be perpendicular, but the resulting array will make
            it look as though they are.
        uvals,vvals
            The nu and nv points along axis U and axis V respectively, 
            that will be used to generate test points for constructing 
            the sample array. Both must be ordered 1D lists or NumPy 
            arrays.
        func
            The 3D function. This function must take a single parameter,
            namely an N*3 array of test points, with N = nu*nv.
    
    RETURNS:
        A 2D NumPy array of sampled function values. This can be used
        directly by Matplotlib's contour() function.
        
    RAISES:
        Can throw NumPy exceptions if the data does not have the
        proper shape.
    
    SEE ALSO:
        makeContourXY()
        makeContourXZ()
        makeContourYZ()
        plotContour()
    """
    p0 = numpy.asarray(p0)
    uaxis = numpy.asarray(uaxis)
    vaxis = numpy.asarray(vaxis)
    uvals = numpy.asarray(uvals)
    vvals = numpy.asarray(vvals)
    nu = uvals.shape[0]
    nv = vvals.shape[0]
    upos = numpy.tile(numpy.c_[uvals,uvals,uvals] * uaxis, (nv,1))
    vpos = numpy.repeat(numpy.c_[vvals,vvals,vvals] * vaxis, nu, axis=0)
    pnts = upos + vpos + p0
    return numpy.reshape(func(pnts), (nv,nu))


def makeContourXY(xvals, yvals, func, dz = 0.0):
    """Prepare data for a contour plot of a 3D function over a planar
    cross section parallel to the XY plane.
    
    PARAMETERS:
        xvals,yvals
            The nx and ny points along the X axis and Y axis respectively, 
            that will be used to generate test points for constructing 
            the sample array. Both must be ordered 1D lists or NumPy 
            arrays.
        func
            The 3D function. This function must take a single parameter,
            namely an N*3 array of test points, with N = nx*ny.
        dz 
            The perpendicular planar increment. Defaults to 0.0, meaning
            the plane crosses the origin.
    
    RETURNS:
        A 2D NumPy array of sampled function values. This can be used
        directly by Matplotlib's contour() function.
    
    RAISES:
        Can throw NumPy exceptions if the data does not have the
        proper shape.
    
    SEE ALSO:
        makeContour()
        makeContourXZ()
        makeContourYZ()
        plotContour()
        matplotlib.contour()
    """
    p0 = numpy.array([0.0, 0.0, dz])
    uaxis = numpy.array([1.0, 0.0, 0.0])
    vaxis = numpy.array([0.0, 1.0, 0.0])
    return makeContour(p0, uaxis, vaxis, xvals, yvals, func)


def makeContourXZ(xvals, zvals, func, dy = 0.0):
    """Prepare data for a contour plot of a 3D function over a planar
    cross section parallel to the XZ plane.
    
    PARAMETERS:
        xvals,zvals
            The nx and nz points along X axis and Z axis respectively, 
            that will be used to generate test points for constructing 
            the sample array. Both must be ordered 1D lists or NumPy 
            arrays.
        func
            The 3D function. This function must take a single parameter,
            namely an N*3 array of test points, with N = nx*nz.
        dy 
            The perpendicular planar increment. Defaults to 0.0, meaning
            the plane crosses the origin.
    
    RETURNS:
        A 2D NumPy array of sampled function values. This can be used
        directly by Matplotlib's contour() function.
    
    RAISES:
        Can throw NumPy exceptions if the data does not have the
        proper shape.
    
    SEE ALSO:
        makeContour()
        makeContourXY()
        makeContourYZ()
        plotContour()
    """
    p0 = numpy.array([0.0, dy, 0.0])
    uaxis = numpy.array([1.0, 0.0, 0.0])
    vaxis = numpy.array([0.0, 0.0, 1.0])
    return makeContour(p0, uaxis, vaxis, xvals, zvals, func)


def makeContourYZ(yvals, zvals, func, dx = 0.0):
    """Prepare data for a contour plot of a 3D function over a planar
    cross section parallel to the YZ plane.
    
    PARAMETERS:
        yvals,zvals
            The ny and nz points along Y axis and Z axis respectively, 
            that will be used to generate test points for constructing 
            the sample array. Both must be ordered 1D lists or NumPy 
            arrays.
        func
            The 3D function. This function must take a single parameter,
            namely an N*3 array of test points, with N = ny*nz.
        dx
            The perpendicular planar increment. Defaults to 0.0, meaning
            the plane crosses the origin.
    
    RETURNS:
        A 2D NumPy array of sampled function values. This can be used
        directly by Matplotlib's contour() function.
    
    RAISES:
        Can throw NumPy exceptions if the data does not have the
        proper shape.
    
    SEE ALSO:
        makeContour()
        makeContourXY()
        makeContourXZ()
        plotContour()
    """
    p0 = numpy.array([dx, 0.0, 0.0])
    uaxis = numpy.array([0.0, 1.0, 0.0])
    vaxis = numpy.array([0.0, 0.0, 1.0])
    return makeContour(p0, uaxis, vaxis, yvals, zvals, func)


def plotContourRB(u,v,data, clip=1.0, d=0.1, lines=True, thin=0.3, thick=3.0):
    """Make a contour plot for a set of 2D signed field data. Negative 
    regions of the field will be colored red, positive regions will be
    colored blue.  
    
    PARAMETERS:
        u,v
           The coordinates along the two axes over which the data has  
           been defined. 
        data
            A 2D array of 'height' data, interpreted over u,v.
        clip
            Limit the colorization to height values in the 
            range [-clip, clip]. Values outside this range are colored
            with on both the positive and negative sides. The default 
            value is 1.0.
        d
            The increment in height value that corresponds to a new 
            level of color; basically the distance in height between
            neighbouring contour lines. Default value is 0.1. 
            (Example: together with clip=1.0, this would mean a total of 
            20 levels of colorization.)
        lines
            A boolean value that controls whether black contour lines
            are drawn on top of the colorized contour map. Defaults 
            to True. 
        thin
             The thickness of the normal contourlines, expressed in points.
             This option is only relevant if option 'lines' is set to True.
             The default value is 0.3, which gives rather fine lines.
        thick
            When contour lines are drawn on top of the colourized 
            contour map, a thicker line is drawn at the 'sea level',
            i.e. locations where the height is 0.0. This option is only
            relevant is option 'lines' is set to True. The default value 
            is 3.0, which gives a rather thick contour line.
    
    RETURNS:
        ---
    """
    # Build colormap
    cdict = { 'red': ( (0.0,1.0,1.0), (0.5,1.0,0.0), (1.0,1.0,1.0) ),
            'green': ( (0.0,1.0,1.0), (0.5,0.0,0.0), (1.0,1.0,1.0) ),
             'blue': ( (0.0,1.0,1.0), (0.5,0.0,1.0), (1.0,1.0,1.0) ) }
    rb = matplotlib.colors.LinearSegmentedColormap('redblue', cdict, 256)
    
    # Prepare data.
    data2 = data.copy()
    data2[data2 < -clip] = -clip
    data2[data2 >  clip] =  clip
    levels = numpy.arange(-clip,clip+d,d)
    
    # Plotting code for gradient
    was_ia = isinteractive()
    ioff()
    contourf(u, v, data2, levels, cmap=rb, hold=True)
    colorbar()
    # Plotting code for level lines
    if lines == True:
        minval = math.floor(data.min() / d) * d
        maxval = math.ceil(data.max() / d) * d
        levels2 = numpy.arange(minval,maxval+d,d)
        contour(u, v, data, levels2, colors='0.0', linewidths=thin)
        levels3 = numpy.array([0.0])
        contour(u, v, data, levels3, colors='0.0', linewidths=thick)
    xlabel('U axis')
    ylabel('V axis')
    axis('image')
    if was_ia == True: 
        ion()
        show()
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
# END
