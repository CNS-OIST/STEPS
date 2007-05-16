# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from steps import IDError


def isValidID(id):
    """Test whether id is a valid identifier for some named STEPS component.
    
    A valid id must be at least 1 character long and must start with an 
    underscore or an alphabetical character (a to z, or A to Z). Following
    characters can be any alphanumerical character or again underscores.
    The id cannot contain spaces.

    Examples of valid id's:
        a, _a, _a_, a000, adasf0, FSDaa9

    Parameters:
        id
            The string that will be tested.

    Returns:
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
    isValidID() returns True or False, checkID() raises a ValueError
    exception if the id is not valid. This makes it useful for use as an
    assertion.

    Parameters:
        id
            The string that will be tested.

    Returns:
        The id itself.

    Raises:
        steps.IDError
            If the id is not valid.
    """
    if not isValidID(id):
        raise IDError, '\'%s\' is not a valid id.' % id
    return id


# END
