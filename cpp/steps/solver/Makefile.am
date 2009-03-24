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

# cpp/steps/solver/Makefile.am

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Include local top source directory to include path.
AM_CPPFLAGS = -I$(top_srcdir)/cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

header_steps_solverdir = $(pkgincludedir)/solver

header_steps_solver_HEADERS = \
    api.hpp \
    doxygen.hpp \
    compdef.hpp \
    diffdef.hpp \
    patchdef.hpp \
    reacdef.hpp \
    specdef.hpp \
    sreacdef.hpp \
    statedef.hpp \
    types.hpp

noinst_LTLIBRARIES = libsteps_solver.la
libsteps_solver_la_SOURCES = \
    api_comp.cpp \
    api_main.cpp \
    api_patch.cpp \
    api_tet.cpp \
    api_tri.cpp \
    compdef.cpp \
    diffdef.cpp \
    patchdef.cpp \
    reacdef.cpp \
    specdef.cpp \
    sreacdef.cpp \
    statedef.cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
