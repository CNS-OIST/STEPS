# cython:language_level=3str
####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

from libcpp.memory cimport shared_ptr
from steps_rng cimport *

# ======================================================================================================================
# Python bindings to namespace steps::rng
# ======================================================================================================================

def _py_rng_create( str rng_name, uint bufsize ):
    """
    Creates and returns a reference to a steps.rng.RNG random number generator object, 
    which is specified by type and pre-allocates a buffer list with size of buffer_size.

    Syntax::
        
        create(type, buffer_size)

    Arguments:
    string type
    int buffer_size

    Return:
    steps.rng.RNG

    """
    return _py_RNG.from_shared_ptr(create(to_std_string(rng_name), bufsize))

def _py_rng_create_mt19937( uint bufsize ):
    """
    Creates and returns a reference to a steps.rng.RNG random number generator object, 
    which is specified by type and pre-allocates a buffer list with size of buffer_size.

    Syntax::
        
        create_mt19937(buffer_size)

    Arguments:
    int buffer_size

    Return:
    steps.rng.RNG

    """
    return _py_RNG.from_shared_ptr(create_mt19937(bufsize))


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_RNG(_py__base):
    "Python wrapper class for RNG"
# ----------------------------------------------------------------------------------------------------------------------
    cdef shared_ptr[RNG] _shared_ptr

    cdef shared_ptr[RNG] ptr(self):
        return self._shared_ptr

    # RNG is abstract
    # def __init__(self, *arg):

    def initialize(self, ulong seed):
        """
        Initialize the random number generator with given seed value.
        
        Syntax::
        
            initialize(seed)
        
        Arguments:
        int seed
        
        Return:
        None
        
        """
        self.ptr().get().initialize(seed)

    def min(self, ):
        """
        Return the minimal inclusive range.
        
        Syntax::
        
        	min()
        
        Arguments:
        None
        
        Return:
        int
        
        """
        return self.get().min()

    def max(self, ):
        """
        Return the maximal inclusive range.
        
        Syntax::
        
        	max()
        
        Arguments:
        None
        
        Return:
        int
        
        """
        return self.get().max()

    def __call__(self, ):
        return deref(self.ptr())()

    def get(self, ):
        """
        Return the next random int in the buffer of the generator.
        
        Syntax::
        
       		get()
        
        Arguments:
        None
        
        Return:
        int
        
        """
        return self.ptr().get().get()

    def getUnfII(self, ):
        """
        Generates a uniform random number on [0,1] real interval.
        
        Syntax::
        
        	getUnfII()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getUnfII()

    def getUnfIE(self, ):
        """
        Generates a uniform random number on [0,1] real interval.
        
        Syntax::
        
        	getUnfIE()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getUnfIE()

    def getUnfEE(self, ):
        """
        Generates a uniform random number on [0,1] real interval.
        
        Syntax::
        
        	getUnfEE()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getUnfEE()

    def getUnfIE53(self, ):
        """
        Generates a uniform random number on [0,1] real interval with 53-bit resolution.
        
        Syntax::
        
        	getUnfIE53()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getUnfIE53()

    def getStdExp(self, ):
        """
        Get a standard exponentially distributed number.
        
        Syntax::
        
        	getStdExp()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getStdExp()

    def getExp(self, double lambda_):
        """
        Get an exponentially distributed number with mean lambda.
        
        Syntax::
        
        	getExp()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getExp(lambda_)

    def getPsn(self, double lambda_):
        """
        Get a Poisson-distributed number with mean lambda.
        
        Syntax::
        
        	getPsn()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getPsn(lambda_)

    def getStdNrm(self, ):
        """
        Get a standard normally distributed random number.
        
        Syntax::
        
        	getStdNrm()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().get().getStdNrm()

    def getBinom(self, uint t, double p):
        """
        Get a binomially distributed number with parameters t and p.
        
        Syntax::
        
        	getBinom()
        
        Arguments:
        int t
        float p
        
        Return:
        int
        
        """
        return self.ptr().get().getBinom(t, p)

    @staticmethod
    cdef _py_RNG from_shared_ptr(shared_ptr[RNG] ptr):
        cdef _py_RNG obj = _py_RNG.__new__(_py_RNG)
        obj._shared_ptr = ptr
        obj._ptr = ptr.get()
        return obj

