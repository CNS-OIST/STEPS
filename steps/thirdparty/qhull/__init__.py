# __init__.py - initialisation script of the Delny package
#
# Copyright 2004 Floris Bruynooghe
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

"""Delaunay triangulation.

This package provides a python interface to the Delaunay triangulation
provided by the qhull software package using the C library libqhull.
It is possible to make a triangulations of any dimension you want.

The Delaunay module imports the core module directy in it's namespace
so you can for example use `Delaunay.Triangulation()' after a `import
Delaunay'.  So see the Delaunay.core documentation for the use of the
module.

Whenever a "sequence object" is mentioned in this package all the
builtin Python sequences are accepted as well as an array from the
Numeric module.
"""

from core import *
from _qhull import delny

__all__ = ["core", "_qhull"]
