###___license_placeholder___###

from steps_rng cimport *

# ======================================================================================================================
# Python bindings to namespace steps::rng
# ======================================================================================================================

def _py_rng_create( str rng_name, unsigned int bufsize ):
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
    return _py_RNG.from_ptr(create(to_std_string(rng_name), bufsize))

def _py_rng_create_mt19937( unsigned int bufsize ):
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
    return _py_RNG.from_ptr(create_mt19937(bufsize))


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_RNG(_py__base):
    "Python wrapper class for RNG"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[RNG] _autodealoc
    cdef RNG *ptr(self):
        return <RNG*> self._ptr

    # RNG is abstract
    # def __init__(self, *arg):

    def initialize(self, unsigned long seed):
        """
        Initialize the random number generator with given seed value.
        
        Syntax::
        
            initialize(seed)
        
        Arguments:
        int seed
        
        Return:
        None
        
        """
        self.ptr().initialize(seed)

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
        return self.ptr().min()

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
        return self.ptr().max()

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
        return self.ptr().get()

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
        return self.ptr().getUnfII()

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
        return self.ptr().getUnfIE()

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
        return self.ptr().getUnfEE()

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
        return self.ptr().getUnfIE53()

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
        return self.ptr().getStdExp()

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
        return self.ptr().getExp(lambda_)

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
        return self.ptr().getPsn(lambda_)

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
        return self.ptr().getStdNrm()

    def getBinom(self, unsigned int t, double p):
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
        return self.ptr().getBinom(t, p)

    @staticmethod
    cdef _py_RNG from_ptr(RNG *ptr):
        cdef _py_RNG obj = _py_RNG.__new__(_py_RNG)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_RNG from_ref(const RNG &ref):
        return _py_RNG.from_ptr(<RNG*>&ref)


