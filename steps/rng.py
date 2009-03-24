# rng.py

"""
This file is the user-interface file for all rng objects in steps.
All objects are directly derived from the corresponding swig objects.
RNG object is owned by Python

"""
####### currently only imports swig file !! ##############

from rng_swig import *

"""
import rng_swig
import _rng_swig

class RNG(rng_swig.RNG) :
    def __init__(self, *args, **kwargs): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _rng_swig.delete_RNG
    __del__ = lambda self : None;
    self.thisown = True

def create_mt19937(*args):
  return rng_swig.create_mt19937(*args)

def create(*args):
  return rng_swig.create(*args)
"""
