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

# cpp/steps/tetexact/Makefile.am

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Include local top source directory to include path.
AM_CPPFLAGS = -I$(top_srcdir)/cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

header_steps_tetexactdir = $(pkgincludedir)/tetexact

header_steps_tetexact_HEADERS = \
    doxygen.hpp \
    comp.hpp \
    diff.hpp \
    kproc.hpp \
    patch.hpp \
    reac.hpp \
    sreac.hpp \
    tet.hpp \
    tri.hpp \
    tetexact.hpp

noinst_LTLIBRARIES = libsteps_tetexact.la
libsteps_tetexact_la_SOURCES = \
    comp.cpp \
    diff.cpp \
    kproc.cpp \
    patch.cpp \
    reac.cpp \
    sreac.cpp \
    tet.cpp \
    tri.cpp \
    tetexact.cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
