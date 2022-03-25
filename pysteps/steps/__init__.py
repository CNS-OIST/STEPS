# -*- coding: utf-8 -*-

####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import atexit
import glob
import importlib.abc
import importlib.machinery
import os.path
import sys
import types
import warnings

__name__      = 'steps'
__longname__  = 'STochastic Engine for Pathway Simulation'
__version__   = '4.0.0'
__author__    = 'STEPS Development Team'
__url__       = 'steps.sourceforge.net'
__license__   = 'GPL2.0'
__binding__   = 'Cython'

#Try importing mpi version if file exists. Avoids catching exception for the wrong reason
if glob.glob(os.path.join(os.path.dirname(__file__),'cysteps_dist.*')):
    from . import cysteps_dist
    from . import cysteps_dist as stepslib
elif glob.glob(os.path.join(os.path.dirname(__file__),'cysteps_mpi.*')):
    from . import cysteps_mpi
    from . import cysteps_mpi as stepslib
elif glob.glob(os.path.join(os.path.dirname(__file__),'cysteps.*')):
    from . import cysteps
    from . import cysteps as stepslib
else:
    raise Exception("STEPS library not found [cysteps*]")

stepslib._py_init()
atexit.register(stepslib._py_finish)

_suppress_greet = False
_quiet = False

def _greet():
    global _suppress_greet
    if not _suppress_greet:
        print("")
        print(__longname__)
        print("Version: ", __version__)
        print("License: ", __license__)
        print("Website: ", __url__)
        print("CXX Binding:", __binding__)

    _suppress_greet = True

###############################
# Custom importing mechanisms #
###############################

class _CustomVirtualLoader(importlib.abc.Loader):
    def create_module(self, spec):
        return types.ModuleType(spec.name)

    def exec_module(self, module):
        pass


class _CustomActualLoader(importlib.abc.Loader):

    def create_module(self, spec):
        return importlib.import_module(spec._actualname)

    # No need for exec since importlib.import_module already takes care of it.
    def exec_module(self, module):
        pass


class _CustomMetaPathFinder(importlib.abc.MetaPathFinder):
    _API_DIRS = [
        'API_1',
        'API_2',
    ]
    _VIRTUAL_IMPORT_PATHS = {
        'steps.interface': 'API_2',
    }
    _DEFAULT_API = 'API_1'

    def __init__(self):
        self._currAPI = _CustomMetaPathFinder._DEFAULT_API
        self._virtualLoader = _CustomVirtualLoader()
        self._actualLoader = _CustomActualLoader()

    def find_spec(self, fullname, path, target=None):
        if fullname in _CustomMetaPathFinder._VIRTUAL_IMPORT_PATHS:
            self._currAPI = _CustomMetaPathFinder._VIRTUAL_IMPORT_PATHS[fullname]
            return importlib.machinery.ModuleSpec(fullname, self._virtualLoader)

        if fullname.startswith(__name__):
            splt = fullname.split('.')
            # If we are trying to load a module from a specific API directory, it means the 
            # default path finders did not succeed loading it, i.e. it does not exist.
            if len(splt) > 1 and splt[1] in _CustomMetaPathFinder._API_DIRS:
                raise ModuleNotFoundError(f'Could not import {fullname}, check that you are '
                                          f'importing modules from the correct API version.')
            newName = '.'.join([splt[0], self._currAPI] + splt[1:])
            # The name of the spec corresponds to the user name (e.g. steps.model) but the actual
            # package to be imported is set in _actualName (e.g. steps.API_1.model). Using directly
            # the actual package path would create issues because the user package path would not
            # be added to sys.modules.
            spec = importlib.machinery.ModuleSpec(fullname, self._actualLoader)
            spec._actualname = newName
            return spec

        return None

sys.meta_path.append(_CustomMetaPathFinder())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
