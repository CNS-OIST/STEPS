# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""This module defines all exceptions that can be raised during execution
of STEPS.

We try to keep the number of distinct exception classes as limited as 
possible, by keeping the scope of each exception class fairly general. 
Instead, we put as much information as possible in the message of the
exception. 
"""


class Error(Exception):
    """Base exception class of the steps package.
    """
    def __init__(self, args):
        self.args = args


class RuntimeError(Exception):
    """Exceptions belonging to this category signal that something is 
    wrong due to external factors, such as user input, or due to 
    the environment (disk full, time out, ...). These should only be
    raised by functions that are directly accessed by users when 
    STEPS is used in a normal fashion.
    """
    def __init__(self, args):
        self.args = args


class ProgramError(Exception):
    """Exceptions in this category typically imply that there
    is a bug. This bug could result from not having checked for
    invalid input arguments earlier. It might also result from a 
    logical error. When an exception of this type is encountered,
    the user should be encouraged to contact the developers.
    """
    def __init__(self, args):
        self.args = args




class IDError(Exception):
    def __init__(self, args):
        self.args = args
class ModelError(Exception):
    def __init__(self, args):
        self.args = args
class GeomError(Exception):
    def __init__(self, args):
        self.args = args






class ArgumentError(RuntimeError):
    """The user specified an argument to some STEPS function that really
    doesn't make sense.
    """
    def __init__(self, args):
        self.args = args


class SolverCoreError(ProgramError):
    """An error occured while loading or working with a solver core module.
    """
    def __init__(self, args):
        self.args = args


# END
