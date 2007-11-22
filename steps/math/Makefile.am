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

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Files that should be cleaned when building a distribution.
distclean-local:
	rm -rf .DS_Store *~ *.bak

pyexec_steps_mathdir = $(pythondir)/steps/math
python_steps_mathdir = $(pythondir)/steps/math
header_steps_mathdir = $(pkgincludedir)/steps/thirdparty/qhull

AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir) \
-I/Library/Frameworks/Python.framework/Versions/2.4/lib/python2.4/site-packages/numpy/core/include

python_steps_math_PYTHON = \
    __init__.py \
    tetrahedron_test.py \
    triangle.py \
    triangle_test.py

pkginclude_HEADERS = \
    accumulate.hpp \
    constants.hpp \
    doxygen.hpp \
    iota.hpp \
    tools.hpp 

pyexec_steps_math_LTLIBRARIES = tetrahedron.la
tetrahedron_la_LDFLAGS = -module -avoid-version $(PYTHON_LDFLAGS)
tetrahedron_la_SOURCES = tetrahedron.cpp
tetrahedron_la_LIBADD = ../rng/libsteps_rng.la
    
# Compiler flags.
AM_CXXFLAGS = \
	-I$(top_srcdir)
    
# END
