# ./steps/sim/shared/Makefile.am
#
# Process this file with automake to produce Makefile.in

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
    reacdef.cpp \
	specdef.cpp \
	statedef.cpp

header_steps_sim_shared_HEADERS = \
	compdef.hpp \
    doxygen.hpp \
    reacdef.hpp \
	specdef.hpp \
	statedef.hpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
