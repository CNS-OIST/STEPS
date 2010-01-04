# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# steps/Makefile.am

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Include local top source directory to include path.
AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir)/cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Process subdirectories.
SUBDIRS = utilities \
          .

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pyexec_stepsdir = $(pythondir)/steps
python_stepsdir = $(pythondir)/steps

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
    geom.py \
    geom_swig.py \
    rng_swig.py \
    rng.py \
    solver_swig.py \
    solver.py


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
