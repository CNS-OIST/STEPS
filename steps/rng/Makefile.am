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

pyexec_steps_rngdir = $(pythondir)/steps/rng
python_steps_rngdir = $(pythondir)/steps/rng
header_steps_rngdir = $(pkgincludedir)/steps/rng

AM_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(PYTHON_CPPFLAGS) -I$(top_srcdir)

SWIG_PYTHON_OPT = @SWIG_PYTHON_OPT@ -I$(top_srcdir)
%_wrap.cpp: %.i
	$(SWIG) $(SWIG_PYTHON_OPT) -o $@ $<

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# LIBRARY: LIBSRNG
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

lib_LTLIBRARIES = libsteps_rng.la

libsteps_rng_la_SOURCES = \
    mt19937.cpp \
    rng.cpp

header_steps_rng_HEADERS = \
    doxygen.hpp \
    mt19937.hpp \
    rng.hpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# EXTENSION: RNG_CORE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

SWIG_SOURCES = rng_core.i

pyexec_steps_rng_LTLIBRARIES = _rng_core.la
_rng_core_la_LDFLAGS = -module -avoid-version $(PYTHON_LDFLAGS)

_rng_core_la_SOURCES = \
    $(SWIG_SOURCES) \
    rng_core_wrap.cpp

_rng_core_la_LIBADD = \
    ./libsteps_rng.la

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PYTHON PACKAGES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

python_steps_rng_PYTHON = \
    __init__.py \
    rng_core.py
    
# END
