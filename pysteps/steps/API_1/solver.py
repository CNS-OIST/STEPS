####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
from __future__ import print_function

from steps import stepslib


# Constants aliases (yep, must be hand coded)
EF_NONE = stepslib._py_API.EF_NONE
EF_DEFAULT = stepslib._py_API.EF_DEFAULT
EF_DV_BDSYS = stepslib._py_API.EF_DV_BDSYS
EF_DV_PETSC  = stepslib._py_API.EF_DV_PETSC


# --------------------------------------------------------------------
# Mixin class which provides helper methods, common to all Solvers
# Method names start with _ to avoid conflicts in multiple inheritance
# --------------------------------------------------------------------
class _Base_Solver(object):
    def _advance_checkpoint_run(self, end_time, cp_interval, prefix, method_name):
        """ Advance the simulation for advance_time,
            automatically checkpoint at each cp_interval.
            Prefix can be added using prefix=<prefix_string>.
            """
        parent = super(self.__class__, self)
        if cp_interval > 0:
            while self.getTime() + cp_interval < end_time:
                self.advance(cp_interval)
                filename = "%s%e.%s_cp" % (prefix, self.getTime(), method_name)
                print("Checkpointing -> ", filename)
                self.checkpoint(filename)

            parent.run(end_time)
            filename = "%s%e.%s_cp" % (prefix, self.getTime(), method_name)
            print("Checkpointing -> ", filename)
            self.checkpoint(filename)
        else:
            parent.run(end_time)

    def _getIndexMapping(self):
        """Get a mapping between compartments/patches/species
           and their indices in the solver."""
        mapping = {"Comp": [], "Patch": []}
        ncomps = self.getNComps()
        for c in range(ncomps):
            cname = self.getCompName(c)
            spec_names = []
            nspces = self.getNCompSpecs(c)
            for s in range(nspces):
                spec_names.append(self.getCompSpecName(c, s))
            mapping["Comp"].append({"Name": cname, "Species": spec_names})

        npatches = self.getNPatches()
        for p in range(npatches):
            pname = self.getPatchName(p)
            spec_names = []
            nspces = self.getNPatchSpecs(p)
            for s in range(nspces):
                spec_names.append(self.getPatchSpecName(p, s))
            mapping["Patch"].append({"Name": pname, "Species": spec_names})
        return mapping


class Wmrk4(stepslib._py_Wmrk4, _Base_Solver):
    """
    Construction::
        
        sim = steps.solver.Wmrk4(model, geom)
        
    Create a non-spatial deterministic solver based on the Runge-Kutta fourth order method.
        
    Arguments:
    steps.model.Model model
    steps.geom.Geom geom
    """
    
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until end_time,
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'wmrk4')

    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for advance_time,
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        end_time = self.getTime() + advance_time
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'wmrk4')
        
    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        return self._getIndexMapping(self)


class Wmdirect(stepslib._py_Wmdirect, _Base_Solver):
    """
    Construction::
    
        sim = steps.solver.Wmdirect(model, geom, rng)
    
    Create a non-spatial stochastic solver based on Gillespie's SSA.
    
    Arguments:
    steps.model.Model model
    steps.geom.Geom geom
    steps.rng.RNG rng
    """
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until <end_time>,
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'wmdirect')
        
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for advance_time,
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        end_time = self.getTime() + advance_time
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'wmdirect')
        
    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        return self._getIndexMapping()


class Wmrssa(stepslib._py_Wmrssa, _Base_Solver):
    """
    Construction::
    
        sim = steps.solver.Wmrssa(model, geom, rng)
    
    Create a non-spatial stochastic solver based on Gillespie's SSA.
    
    Arguments:
    steps.model.Model model
    steps.geom.Geom geom
    steps.rng.RNG rng
    """
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until <end_time>,
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'wmrssa')
        
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for advance_time,
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        end_time = self.getTime() + advance_time
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'wmrssa')
        
    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        return self._getIndexMapping()
        
        
class Tetexact(stepslib._py_Tetexact, _Base_Solver):
    """
    Construction::
    
        sim = steps.solver.Tetexact(model, geom, rng, calcMembPot = 0)
    
    Create a spatial stochastic solver based on Gillespie's SSA, extended with diffusion across elements in a tetrahedral mesh.
    If voltage is to be simulated, argument calcMemPot=1 will set to the default solver. calcMembPot=0 means voltage will not be simulated. 
    
    Arguments:
    steps.model.Model model
    steps.geom.Geom geom
    steps.rng.RNG rng
    int calcMemPot (default=0)
    
    """
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until <end_time>,
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'tetexact')
        
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for <advance_time>,
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        end_time = self.getTime() + advance_time
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'tetexact')
        
    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        return self._getIndexMapping()


class TetODE(stepslib._py_TetODE, _Base_Solver):
    """
    Construction::
    
        sim = steps.solver.TetODE(model, geom, rng=None, calcMembPot = 0)
    
    Create a spatial determinstic solver based on the CVODE library.
    If voltage is to be simulated, argument calcMemPot=1 will set to the default solver. calcMembPot=0 means voltage will not be simulated. 
    
    Arguments:
    steps.model.Model model
    steps.geom.Geom geom
    steps.rng.RNG rng (default=None)
    int calcMemPot (default=0)
    
    """
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until <end_time>,
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'tetode')
        
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for <advance_time>,
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        end_time = self.getTime() + advance_time
        self._advance_checkpoint_run(end_time, cp_interval, prefix, 'tetode')
    
    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        return self._getIndexMapping()
