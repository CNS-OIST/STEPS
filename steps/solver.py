# solver.py

"""

This file is the user-interface file for all solver objects.
All objects are directly derived from the corresponding swig objects.
All objects are owned by Python.

"""

import solver_swig
import _solver_swig

class Wmrk4(solver_swig.Wmrk4) :
    def __init__(self, *args): 
        """__init__(self, steps::model::Model m, steps::wm::Geom g, steps::rng::RNG r) -> Wmrk4"""
        this = _solver_swig.new_Wmrk4(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True

class Wmdirect(solver_swig.Wmdirect) :
    def __init__(self, *args): 
        """__init__(self, steps::model::Model m, steps::wm::Geom g, steps::rng::RNG r) -> Wmdirect"""
        this = _solver_swig.new_Wmdirect(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True

class Tetexact(solver_swig.Tetexact) :
    def __init__(self, *args): 
        """__init__(self, steps::model::Model m, steps::wm::Geom g, steps::rng::RNG r) -> Tetexact"""
        this = _solver_swig.new_Tetexact(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = True