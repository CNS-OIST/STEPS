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

header_steps_sim_shareddir = $(pkgincludedir)/steps/sim/shared

AM_CPPFLAGS = -I$(top_srcdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# LIBRARY: LIBSWMDIRECT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

noinst_LTLIBRARIES = libsimshared.la

libsimshared_la_SOURCES = \
	compdef.cpp \
    diffdef.cpp \
    reacdef.cpp \
	specdef.cpp \
	statedef.cpp

header_steps_sim_shared_HEADERS = \
	compdef.hpp \
    diffdef.hpp \
    doxygen.hpp \
    reacdef.hpp \
	specdef.hpp \
	statedef.hpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
