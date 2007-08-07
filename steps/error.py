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


"""This module defines all exceptions that can be raised within STEPS code.

We try to keep the number of distinct exception classes as limited as 
possible, by keeping the scope of each exception class fairly general. 
Instead, we put as much information as possible in the message of the
exception. 
"""


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Error(Exception):
    """Base exception class of the steps package.
    """
    def __init__(self, args):
        self.args = args


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class RuntimeError(Exception):
    """Exceptions belonging to this category signal that something is 
    wrong due to external factors, such as user input, or due to 
    the environment (disk full, time out, ...). These should only be
    raised by functions that are directly accessed by users when 
    STEPS is used in a normal fashion.
    """
    def __init__(self, args):
        self.args = args


class ArgumentError(RuntimeError):
    """The user specified an argument to some STEPS function that somehow
    doesn't make sense. 
    
    This is one of the most commonly raised exceptions.
    """
    def __init__(self, args):
        self.args = args


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class ProgramError(Exception):
    """Exceptions in this category typically imply that there
    is a bug. This bug could result from not having checked for
    invalid input arguments earlier. It might also result from a 
    logical error. When an exception of this type is encountered,
    the user should be encouraged to contact the developers.
    """
    def __init__(self, args):
        self.args = args


class SolverCoreError(ProgramError):
    """An error occured while loading or working with a solver core module.
    """
    def __init__(self, args):
        self.args = args


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
