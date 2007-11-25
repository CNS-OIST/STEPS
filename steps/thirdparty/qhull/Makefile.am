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

# Set $MAKE if necessary.
@SET_MAKE@

# VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

# Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak

pyexec_steps_thirdparty_qhulldir = $(pythondir)/steps/thirdparty/qhull
python_steps_thirdparty_qhulldir = $(pythondir)/steps/thirdparty/qhull
header_steps_thirdparty_qhulldir = $(pkgincludedir)/steps/thirdparty/qhull

AM_CPPFLAGS = $(PYTHON_CPPFLAGS) $(NUMPY_CPPFLAGS) -I$(top_srcdir)

pyexec_steps_thirdparty_qhull_LTLIBRARIES = _qhull.la
_qhull_la_LDFLAGS = -module -avoid-version $(PYTHON_LDFLAGS)

_qhull_la_SOURCES = \
    _qhullmodule.c \
    geom.c \
    geom2.c \
    global.c \
    io.c \
    mem.c \
    merge.c \
    poly.c \
    poly2.c \
    qhull.c \
    qset.c \
    stat.c \
    user.c

header_steps_thirdparty_qhull_HEADERS = \
    geom.h \
    io.h \
    mem.h \
    merge.h \
    poly.h \
    qhull.h \
    qhull_a.h \
    qset.h \
    stat.h \
    user.h

python_steps_thirdparty_qhull_PYTHON = \
    __init__.py \
    core.py
    
# END
