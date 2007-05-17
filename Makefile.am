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
    doc \
    extras \
    steps \
    thirdparty \
    .

# Things that should go in the distribution.
EXTRA_DIST = \
    setup.py

# Headers that are generated during the build but should be installed.
nobase_nodist_include_HEADERS = \
    steps/config.h

# Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak

# END
