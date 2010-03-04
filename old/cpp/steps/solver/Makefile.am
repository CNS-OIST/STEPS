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

# cpp/steps/solver/Makefile.am

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Include local top source directory to include path.
AM_CPPFLAGS = -I$(top_srcdir)/cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

header_steps_solverdir = $(pkgincludedir)/solver

header_steps_solver_HEADERS = \
    api.hpp \
    doxygen.hpp \
    compdef.hpp \
    diffdef.hpp \
    patchdef.hpp \
    reacdef.hpp \
    specdef.hpp \
    sreacdef.hpp \
    statedef.hpp \
    types.hpp

noinst_LTLIBRARIES = libsteps_solver.la
libsteps_solver_la_SOURCES = \
    api_comp.cpp \
    api_main.cpp \
    api_patch.cpp \
    api_tet.cpp \
    api_tri.cpp \
    compdef.cpp \
    diffdef.cpp \
    patchdef.cpp \
    reacdef.cpp \
    specdef.cpp \
    sreacdef.cpp \
    statedef.cpp

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
