# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
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
