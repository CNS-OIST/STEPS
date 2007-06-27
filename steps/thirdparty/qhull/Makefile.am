# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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

AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir) \
-I/Library/Frameworks/Python.framework/Versions/2.4/lib/python2.4/site-packages/numpy/core/include

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
