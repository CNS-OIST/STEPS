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

# Process subdirectories.
SUBDIRS	= \
    model \
    sim \
    .

# Python modules in local packages.
steps_testsdir = $(pythondir)/steps/tests
steps_tests_PYTHON = \
    __init__.py \
    test.py \
    test_misc.py \
    tools.py

# Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
# END
