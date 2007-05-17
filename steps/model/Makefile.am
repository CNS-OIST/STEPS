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
python_steps_modeldir = $(pythondir)/steps/model
python_steps_model_PYTHON = \
    __init__.py \
    model.py \
    reaction.py \
    species.py \
    volsys.py

# Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
# END
