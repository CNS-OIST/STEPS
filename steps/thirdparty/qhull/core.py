# core.py - the main user interface to the Delny package
#
# Copyright 2004-2006 Floris Bruynooghe
#
# This file is part of Delny.
#
# Delny is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Delny is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Delny; if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
#
#
# Authors: Floris Bruynooghe (flub)

"""The main user interface to the Delny package

Most users will want to use the Triangulate class to access the
routines provided in the package.
"""

import numpy
import _qhull


class Triangulation:
    """Represents a Delaunay triangulation of a set of points

    All sequences used in this class can be either Python builtin
    sequences or Numeric sequences.  Results are returned as Python
    builtin sequences however.
    """
    # Data attributes:
    #   self.neighbours
    #   self.facets
    #   self.indices
    def __init__(self, inputset, dim=None):
        """Creates the Delaunay triangulation

        The `set' can be a 2 dimensional sequence with the last
        dimension being the coordinates of every point.

        Alternatively the set is a sequence objet of any dimension.
        In this case the `dim' argument should be given and indicate
        the number of dimensions.  In this case the sequence will be
        flattened (made 1-dimensional) and the first point will be
        represented by the `dim' first elements of it and so on.
        """
        if(dim != None):
            self.set = Numeric.array(inputset)
            self.set.shape = (-1, dim)
        else:
            self.set = inputset
        self.neighbours, self.facets, self.indices = _qhull.delny(self.set)

    def get_set(self):
        """Returns the set as it is being used by this class

        This could be any sequence object, however if the `dim'
        argument was not passed along to the constructor of the object
        you can be sure this is the same sequence object as you passed
        to the constructor.
        """
        return self.set

    def get_neighbours(self):
        """Returns the neighbours of each point in the set

        This is a dictionnary with the points of the sets as key and a
        list of it's nearest neighbours as value.  Every neighbour is
        a tuple with <dimension> floats.
        """
        return self.neighbours

    def get_elements(self):
        """Returns the elements of the Delaunay triangulation

        This is a list of elements where every element is a list of
        nodes and every node is a tuple of <dimesnion> floats.  An
        element is a triangle in 2D and a tetraheron in 3D.
        """
        return self.facets

    def get_elements_indices(self):
        """Returns the elements of the Delaunay triangulation

        This is a list of elements where every element is a list of
        node indices corresponding to the point index given in the inputset.
        An element is a triangle in 2D and a tetraheron in 3D.
        """
        return self.indices


    def update_set(self, newset):
        """Recalculate the neighbours with a new input set

        This has the same effect as creating a new instance but
        without doing so.
        """
        #FIXME: should this be renamed to set_set()?
        #FIXME: should this also take the `dim' argument?
        self.__init__(newset)
