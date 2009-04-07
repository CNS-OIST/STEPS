# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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

# cpp/steps/geom/Makefile.am

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Include local top source directory to include path.
AM_CPPFLAGS = -I$(top_srcdir)/cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

header_steps_geomdir = $(pkgincludedir)/geom

header_steps_geom_HEADERS = \
    doxygen.hpp \
    comp.hpp \
    geom.hpp \
    patch.hpp \
    tet.hpp \
    tetmesh_rw.hpp \
    tetmesh.hpp \
    tmcomp.hpp \
    tmpatch.hpp \
    tri.hpp

noinst_LTLIBRARIES = libsteps_geom.la
libsteps_geom_la_SOURCES = \
    comp.cpp \
    geom.cpp \
    patch.cpp \
    tet.cpp \
    tetmesh_rw.cpp \
    tetmesh.cpp \
    tmcomp.cpp \
    tmpatch.cpp \
    tri.cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
