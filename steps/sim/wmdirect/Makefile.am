# ./steps/sim/wmdirect/Makefile.am
#
# Process this file with automake to produce Makefile.in

# $Id$

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@SET_MAKE@

srcdir = @srcdir@
VPATH = @srcdir@

distclean-local:
	-rm -rf .DS_Store *~ *.bak

header_steps_sim_wmdirectdir = $(pkgincludedir)/steps/sim/wmdirect

AM_CPPFLAGS = -I$(top_srcdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# LIBRARY: LIBSWMDIRECT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

lib_LTLIBRARIES = libsteps_wmdirect.la

libsteps_wmdirect_la_SOURCES = \
	func_core.cpp \
    func_ssa.cpp \
    state.cpp

libsteps_wmdirect_la_LIBADD = \
    ../shared/libsimshared.la \
	../../rng/libsteps_rng.la

header_steps_sim_wmdirect_HEADERS = \
	doxygen.hpp \
	state.hpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
