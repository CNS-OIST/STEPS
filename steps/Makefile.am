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

# steps/Makefile.am

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
#SUBDIRS = 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pyexec_stepsdir = $(pythondir)/steps_0_5_0
python_stepsdir = $(pythondir)/steps_0_5_0

pyexec_steps_LTLIBRARIES = \
    _model_swig.la \
    _geom_swig.la \
    _rng_swig.la \
    _solver_swig.la

_model_swig_la_LDFLAGS = -module -avoid-version -no-undefined $(PYTHON_LDFLAGS) 
_model_swig_la_SOURCES = model_wrap.cpp
_model_swig_la_LIBADD = $(top_srcdir)/cpp/steps/libsteps.la

#_error_la_LDLFAGS = -module -avoid-version -no-undefined $(PYTHON_LDFLAGS)
#_error_la_SOURCES = error_wrap.cpp
#_error_la_LIBADD = $(top_srcdir)/cpp/steps/libsteps.la

_geom_swig_la_LDFLAGS = -module -avoid-version -no-undefined $(PYTHON_LDFLAGS) 
_geom_swig_la_SOURCES = geom_wrap.cpp
_geom_swig_la_LIBADD = $(top_srcdir)/cpp/steps/libsteps.la

_rng_swig_la_LDFLAGS = -module -avoid-version -no-undefined $(PYTHON_LDFLAGS)
_rng_swig_la_SOURCES = rng_wrap.cpp
_rng_swig_la_LIBADD = $(top_srcdir)/cpp/steps/libsteps.la

_solver_swig_la_LDFLAGS = -module -avoid-version -no-undefined $(PYTHON_LDFLAGS)
_solver_swig_la_SOURCES = solver_wrap.cpp
_solver_swig_la_LIBADD = $(top_srcdir)/cpp/steps/libsteps.la

python_steps_PYTHON = \
    __init__.py \
    model.py \
    model_swig.py \
	error_swig.py \
    geom.py \
    geom_swig.py \
    rng_swig.py \
    rng.py \
    solver_swig.py \
    solver.py \
    tetgen.py

# TODO: add model.i (or *.i) to package, even though it's not realy used.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
