# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
# 
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""
"""


from rng_core import *


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def create(type, bufsize):
    """
    """
    from steps.error import Error
    if type == 'mt':
        return rng_core.create_mt19937(bufsize)
    elif type == 'mt19937':
        return rng_core.create_mt19937(bufsize)
    else:
        raise Error, '\'%s\' is not a valid random number generator' % ( type )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
