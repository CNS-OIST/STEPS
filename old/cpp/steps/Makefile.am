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

# cpp/steps/Makefile.am

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Include local top source directory to include path.
AM_CPPFLAGS = -I$(top_srcdir)/cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Process subdirectories.
SUBDIRS = \
    compat \
    console \
    geom \
    math \
    model \
    rng \
    solver \
    tools \
    wmrk4 \
    tetexact \
    wmdirect \
    .

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

header_stepsdir = $(pkgincludedir)

# Headers that should be included in the distribution and installed.
header_steps_HEADERS = \
    common.h \
    error.hpp

# STEPS library.
lib_LTLIBRARIES = libsteps.la
libsteps_la_SOURCES = \
    steps.cpp \
    error.cpp
libsteps_la_LIBADD = \
    compat/libsteps_compat.la \
    console/libsteps_console.la \
    math/libsteps_math.la \
    model/libsteps_model.la \
    geom/libsteps_geom.la \
    rng/libsteps_rng.la \
    solver/libsteps_solver.la \
    wmrk4/libsteps_wmrk4.la \
    tetexact/libsteps_tetexact.la \
    wmdirect/libsteps_wmdirect.la

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
