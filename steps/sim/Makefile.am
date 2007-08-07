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

pyexec_steps_simdir = $(pythondir)/steps/sim
python_steps_simdir = $(pythondir)/steps/sim
header_steps_simdir = $(pkgincludedir)/steps/sim

header_steps_sim_HEADERS = \
    doxygen.hpp

AM_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(PYTHON_CPPFLAGS) -I$(top_srcdir)

SWIG_PYTHON_OPT = @SWIG_PYTHON_OPT@ -I$(top_srcdir) 
%_wrap.cpp: %.i
	$(SWIG) $(SWIG_PYTHON_OPT) -o $@ $<
SWIG_SOURCES = tetexact_core.i wmdirect_core.i


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# EXTENSIONS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pyexec_steps_sim_LTLIBRARIES = \
    _tetexact_core.la \
    _wmdirect_core.la

# TETEXACT_CORE
_tetexact_core_la_LDFLAGS = -module -avoid-version $(PYTHON_LDFLAGS)
_tetexact_core_la_SOURCES = \
    $(SWIG_SOURCES) \
    tetexact_core_wrap.cpp
_tetexact_core_la_LIBADD = \
    ./tetexact/libsteps_tetexact.la

# WMDIRECT_CORE
_wmdirect_core_la_LDFLAGS = -module -avoid-version $(PYTHON_LDFLAGS)
_wmdirect_core_la_SOURCES = \
    $(SWIG_SOURCES) \
    wmdirect_core_wrap.cpp
_wmdirect_core_la_LIBADD = \
    ./wmdirect/libsteps_wmdirect.la  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PYTHON PACKAGES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

python_steps_sim_PYTHON = \
	__init__.py \
	controller.py \
	tetexact.py \
	tetexact_core.py \
	wmdirect.py \
	wmdirect_core.py
	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

SUBDIRS = \
    shared \
    swiginf \
    tetexact \
	wmdirect \
	.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
