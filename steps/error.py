# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Error(Exception):
    """Base exception class of the steps package."""
    def __init__(self, args):
        self.args = args


class IDError(Error):
    """Invalid identifier for a component."""
    def __init__(self, args):
        self.args = args


class ModelError(Error): 
    """A operation on one of the model components failed."""
    def __init__(self, args):
        self.args = args


class GeomError(Error): 
    """A operation on one of the geometry components failed."""
    def __init__(self, args):
        self.args = args


# END
