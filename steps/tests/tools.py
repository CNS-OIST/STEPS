# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
import sys

def add_steps_path():
    """Make sure all unit tests can reference to steps modules, regardless
    of where they are located."""
    curdir = os.path.dirname(os.path.abspath(__file__)) 
    stepsdir = os.path.normpath(os.path.join(curdir, '..', '..'))
    if not os.path.exists(stepsdir):
        sys.exit(-1)
    if stepsdir not in sys.path:
        sys.path.insert(0, stepsdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
