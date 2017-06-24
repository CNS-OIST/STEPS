# -*- coding: utf-8 -*-

####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# This file is the user-interface file for all solver objects.
# All objects are directly derived from the corresponding swig objects.
# All objects are owned by Python.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
"""
Implementation of simulation solvers. 

Each solver is a partial or full implementation of the STEPS solver API.
At the moment STEPS implements three different solvers. 

The steps.solver.Wmrk4 class implements a well-mixed, deterministic solver 
based on the Runge–Kutta method. 

The steps.solver.Wmdirect class implements a stochastic, well-mixed solver 
based on Gillespie's Direct SSA Method. 

The steps.solver.Tetexact class implements a stochastic reaction-diffusion 
solver, based on Gillespie's Direct SSA Method extended for diffusive fluxes 
between tetrahedral elements in complex geometries.

"""
try:
    from . import steps_swig_numpy as steps_swig
    import _steps_swig_numpy as _steps_swig
except:
    from . import steps_swig
    import _steps_swig

import steps
steps._greet()

import cPickle

# Better way to do this?

EF_NONE = steps_swig.EF_NONE
EF_DEFAULT = steps_swig.EF_DEFAULT
EF_DV_BDSYS = steps_swig.EF_DV_BDSYS
EF_DV_SLUSYS = steps_swig.EF_DV_SLUSYS
EF_DV_PETSC = steps_swig.EF_DV_PETSC

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Well-mixed RK4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Wmrk4(steps_swig.Wmrk4) :
    def __init__(self, model, geom, rng=None): 
        """
        Construction::
        
            sim = steps.solver.Wmrk4(model, geom, rng)
            
        Create a well-mixed RK4 simulation solver.
            
        Arguments: 
            * steps.model.Model model
            * steps.geom.Geom geom
            * steps.rng.RNG rng
        """
        this = _steps_swig.new_Wmrk4(model, geom, rng)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 1
        self.model = model
        self.geom = geom
        
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until end_time, 
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.wmrk4_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.wmrk4_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)
        
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for advance_time, 
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        
        end_time = _steps_swig.API_getTime(self) + advance_time
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.wmrk4_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.wmrk4_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)

    def getIndexMapping(self):
        """
        Get a mapping between compartments/patches/species
        and their indices in the solver.
        """
        mapping = {"Comp":[], "Patch":[]}
        ncomps = _steps_swig.API_getNComps(self)
        for c in range(ncomps):
            cname = _steps_swig.API_getCompName(self, c)
            spec_names = []
            nspces = _steps_swig.API_getNCompSpecs(self, c)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getCompSpecName(self, c, s))
            mapping["Comp"].append({"Name":cname, "Species":spec_names})
        npatches = _steps_swig.API_getNPatches(self)
        for p in range(npatches):
            pname = _steps_swig.API_getPatchName(self, p)
            spec_names = []
            nspces = _steps_swig.API_getNPatchSpecs(self, p)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getPatchSpecName(self, p, s))
            mapping["Patch"].append({"Name":cname, "Species":spec_names})
        return mapping


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Well-mixed Direct SSA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Wmdirect(steps_swig.Wmdirect) :
    def __init__(self, model, geom, rng): 
        """
        Construction::
        
            sim = steps.solver.Wmdirect(model, geom, rng)
            
        Create a well-mixed Direct SSA simulation solver.
            
        Arguments: 
            * steps.model.Model model
            * steps.geom.Geom geom
            * steps.rng.RNG rng
        """
        this = _steps_swig.new_Wmdirect(model, geom, rng)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 1
        self.model = model
        self.geom = geom
        
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until <end_time>, 
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.wmdirect_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.wmdirect_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)
        
    def advance(self, advance_time, cp_interval = 0.0):
        """
        Advance the simulation for advance_time, 
        automatically checkpoint at each cp_interval.
        Prefix can be added using prefix=<prefix_string>.
        """
        
        end_time = _steps_swig.API_getTime(self) + advance_time
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.wmdirect_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.wmdirect_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)
            
    def getIndexMapping(self):
        """
            Get a mapping between compartments/patches/species
            and their indices in the solver.
            """
        mapping = {"Comp":[], "Patch":[]}
        ncomps = _steps_swig.API_getNComps(self)
        for c in range(ncomps):
            cname = _steps_swig.API_getCompName(self, c)
            spec_names = []
            nspces = _steps_swig.API_getNCompSpecs(self, c)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getCompSpecName(self, c, s))
            mapping["Comp"].append({"Name":cname, "Species":spec_names})
        npatches = _steps_swig.API_getNPatches(self)
        for p in range(npatches):
            pname = _steps_swig.API_getPatchName(self, p)
            spec_names = []
            nspces = _steps_swig.API_getNPatchSpecs(self, p)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getPatchSpecName(self, p, s))
            mapping["Patch"].append({"Name":cname, "Species":spec_names})
        return mapping

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Tetrahedral Direct SSA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #        
class Tetexact(steps_swig.Tetexact) :
    def __init__(self, model, geom, rng, calcMembPot = False):
        """
            Construction::
            
            sim = steps.solver.Tetexact(model, geom, rng, calcMembPot = False)
            
            Create a Tetexact SSA simulation solver.
            
            Arguments:
            * steps.model.Model model
            * steps.geom.Geom geom
            * steps.rng.RNG rng
            # int calcMembPot
            
            """
        this = _steps_swig.new_Tetexact(model, geom, rng, calcMembPot)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 1
        self.model = model
        self.geom = geom

        
    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
        Run the simulation until <end_time>, 
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.tetexact_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.tetexact_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)
        
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
        Advance the simulation for <advance_time>, 
        automatically checkpoint at each <cp_interval>.
        Prefix can be added using prefix=<prefix_string>.
        """
        
        end_time = _steps_swig.API_getTime(self) + advance_time
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.tetexact_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.tetexact_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)

    def getIndexMapping(self):
        """
            Get a mapping between compartments/patches/species
            and their indices in the solver.
            """
        mapping = {"Comp":[], "Patch":[]}
        ncomps = _steps_swig.API_getNComps(self)
        for c in range(ncomps):
            cname = _steps_swig.API_getCompName(self, c)
            spec_names = []
            nspces = _steps_swig.API_getNCompSpecs(self, c)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getCompSpecName(self, c, s))
            mapping["Comp"].append({"Name":cname, "Species":spec_names})
        npatches = _steps_swig.API_getNPatches(self)
        for p in range(npatches):
            pname = _steps_swig.API_getPatchName(self, p)
            spec_names = []
            nspces = _steps_swig.API_getNPatchSpecs(self, p)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getPatchSpecName(self, p, s))
            mapping["Patch"].append({"Name":cname, "Species":spec_names})
        return mapping


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Tetrahedral-based ODE solver
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #        
class TetODE(steps_swig.TetODE) :  
    def __init__(self, model, geom, rng=None, calcMembPot = False): 
        """
            Construction::
            
            sim = steps.solver.TetODE(model, geom, rng=0)
            
            Create a TetODE determinstic diffusion solver.
            
            Arguments: 
            * steps.model.Model model
            * steps.geom.Geom geom
            * steps.rng.RNG rng
            # int calcMemPot
            """
        this = _steps_swig.new_TetODE(model, geom, rng, calcMembPot)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 1
        self.model = model
        self.geom = geom

    def run(self, end_time, cp_interval = 0.0, prefix = ""):
        """
            Run the simulation until <end_time>, 
            automatically checkpoint at each <cp_interval>.
            Prefix can be added using prefix=<prefix_string>.
            """
        
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.tetode_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.tetode_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)
    
    def advance(self, advance_time, cp_interval = 0.0, prefix = ""):
        """
            Advance the simulation for <advance_time>, 
            automatically checkpoint at each <cp_interval>.
            Prefix can be added using prefix=<prefix_string>.
            """
        
        end_time = _steps_swig.API_getTime(self) + advance_time
        if cp_interval > 0:
            while _steps_swig.API_getTime(self) + cp_interval < end_time:
                _steps_swig.API_advance(self, cp_interval)
                filename = "%s%e.tetode_cp" % (prefix, _steps_swig.API_getTime(self))
                print "Checkpointing -> ", filename
                _steps_swig.API_checkpoint(self, filename)
            _steps_swig.API_run(self, end_time)
            filename = "%s%e.tetode_cp" % (prefix, _steps_swig.API_getTime(self))
            print "Checkpointing -> ", filename
            _steps_swig.API_checkpoint(self, filename)
        else:
            _steps_swig.API_run(self, end_time)

    def getIndexMapping(self):
        """
            Get a mapping between compartments/patches/species
            and their indices in the solver.
            """
        mapping = {"Comp":[], "Patch":[]}
        ncomps = _steps_swig.API_getNComps(self)
        for c in range(ncomps):
            cname = _steps_swig.API_getCompName(self, c)
            spec_names = []
            nspces = _steps_swig.API_getNCompSpecs(self, c)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getCompSpecName(self, c, s))
            mapping["Comp"].append({"Name":cname, "Species":spec_names})
        npatches = _steps_swig.API_getNPatches(self)
        for p in range(npatches):
            pname = _steps_swig.API_getPatchName(self, p)
            spec_names = []
            nspces = _steps_swig.API_getNPatchSpecs(self, p)
            for s in range(nspces):
                spec_names.append(_steps_swig.API_getPatchSpecName(self, p, s))
            mapping["Patch"].append({"Name":cname, "Species":spec_names})
        return mapping
