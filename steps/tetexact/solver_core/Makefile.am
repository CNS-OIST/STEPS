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

@SET_MAKE@

srcdir = @srcdir@
VPATH = @srcdir@

distclean-local:
	-rm -rf .DS_Store *~ *.bak

header_steps_tetexact_solver_coredir = \
    $(pkgincludedir)/steps/tetexact/solver_core

AM_CPPFLAGS = -I$(top_srcdir)

lib_LTLIBRARIES = libsteps_tetexact.la

libsteps_tetexact_la_SOURCES = \
    func_core.cpp \
    func_ssa.cpp \
    state.cpp

libsteps_tetexact_la_LIBADD = \
    ../../sim/shared/libsimshared.la \
    ../../rng/libsteps_rng.la

header_steps_tetexact_solver_core_HEADERS = \
    doxygen.hpp \
    state.hpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
