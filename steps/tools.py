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


"""A collection of functions that are useful in widely different parts 
of STEPS.
"""


import steps.error as serr


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def isValidID(id):
    """Test whether id is a valid identifier for some named STEPS component.
    
    A valid id must be at least 1 character long and must start with an 
    underscore or an alphabetical character (a to z, or A to Z). Following
    characters can be any alphanumerical character or again underscores.
    The id cannot contain spaces.

    Examples of valid id's:
        a, _a, _a_, a000, adasf0, FSDaa9

    PARAMETERS:
        id
            The string that will be tested.

    RETURNS:
        True
            If the id is valid.
        False
            Otherwise.
    """
    from re import match
    if match(r'[a-zA-Z_][a-zA-Z0-9_]*$', id) == None: return False
    return True


def checkID(id):
    """Test whether id is a valid identifier for some named STEPS component.

    This function calls steps.model.isValidID() for the test. But whereas
    isValidID() returns True or False, checkID() raises a 
    steps.error.ArgumentError exception if the id is not valid. This makes 
    it useful for use as an assertion.

    PARAMETERS:
        id
            The string that will be tested.

    RETURNS:
        The id itself.

    RAISES:
        steps.error.ArgumentError
            If the id is not valid.
    """
    if not isValidID(id):
        raise serr.ArgumentError, '\'%s\' is not a valid id.' % id
    return id


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# END
