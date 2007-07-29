# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@SET_MAKE@

srcdir = @srcdir@
VPATH = @srcdir@

distclean-local:
	-rm -rf .DS_Store *~ *.bak

header_steps_sim_tetexactdir = $(pkgincludedir)/steps/sim/tetexact

AM_CPPFLAGS = -I$(top_srcdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# LIBRARY: LIBSTEPS_TETEXACT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

lib_LTLIBRARIES = libsteps_tetexact.la

libsteps_tetexact_la_SOURCES = \
	func_core.cpp \
    func_ssa.cpp \
    state.cpp

libsteps_tetexact_la_LIBADD = \
    ../shared/libsimshared.la \
	../../rng/libsteps_rng.la

header_steps_sim_tetexact_HEADERS = \
	doxygen.hpp \
	state.hpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
