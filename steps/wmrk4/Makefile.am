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
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@SET_MAKE@

srcdir = @srcdir@
VPATH = @srcdir@

distclean-local:
	-rm -rf .DS_Store *~ *.bak

pyexec_steps_wmrk4dir = $(pythondir)/steps/wmrk4
python_steps_wmrk4dir = $(pythondir)/steps/wmrk4
header_steps_wmrk4dir = $(pkgincludedir)/steps/wmrk4

AM_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(PYTHON_CPPFLAGS) -I$(top_srcdir)

SWIG_PYTHON_OPT = @SWIG_PYTHON_OPT@ -I$(top_srcdir) 
%_wrap.cpp: %.i
	$(SWIG) $(SWIG_PYTHON_OPT) -o $@ $<
SWIG_SOURCES = solver_core.i

pyexec_steps_wmrk4_LTLIBRARIES = \
    _solver_core.la

_solver_core_la_LDFLAGS = -module -avoid-version -no-undefined $(PYTHON_LDFLAGS)
_solver_core_la_SOURCES = \
    $(SWIG_SOURCES) \
    solver_core_wrap.cpp
_solver_core_la_LIBADD = \
    ./solver_core/libsteps_wmrk4.la

python_steps_wmrk4_PYTHON = \
    __init__.py \
    solver.py \
    solver_core.py

SUBDIRS = \
    solver_core \
    .

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
