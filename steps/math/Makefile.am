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

# Python modules in local packages.
python_steps_mathdir = $(pythondir)/steps/math
python_steps_math_PYTHON = \
    __init__.py \
    tetrahedron.py \
    tetrahedron_test.py \
    triangle.py \
    triangle_test.py

# Headers that should be included in the distribution and installed.
pkginclude_HEADERS = \
    accumulate.hpp \
    constants.hpp \
    doxygen.hpp \
    iota.hpp \
    linalg.hpp \
    tools.hpp 

# Compiler flags.
AM_CXXFLAGS = \
	-I$(top_srcdir)

# Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
# END
