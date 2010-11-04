# -*- coding: utf-8 -*-

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev$
#  Last Changed Date: $Date$
#  Last Changed By:   $Author$

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
from . import steps_swig
import _steps_swig
import cPickle

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Well-mixed RK4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Wmrk4(steps_swig.Wmrk4) :
    def __init__(self, model, geom, rng): 
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
        self.thisown = True
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
        self.thisown = True
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
            
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Tetrahedral Direct SSA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #        
class Tetexact(steps_swig.Tetexact) :  
    def __init__(self, model, geom, rng): 
        """
        Construction::
        
            sim = steps.solver.Tetexact(model, geom, rng)
            
        Create a Tetexact SSA simulation solver.
            
        Arguments: 
            * steps.model.Model model
            * steps.geom.Geom geom
            * steps.rng.RNG rng
        """
        print "make solver"
        this = _steps_swig.new_Tetexact(model, geom, rng)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True
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
