###___license_placeholder___###

from steps_rng cimport *

# ======================================================================================================================
# Python bindings to namespace steps::rng
# ======================================================================================================================

def _py_rng_create( std.string rng_name, unsigned int bufsize ):
    return _py_RNG.from_ptr(create(rng_name, bufsize))

def _py_rng_create_mt19937( unsigned int bufsize ):
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
        self.ptr().initialize(seed)

    def min(self, ):
        return self.ptr().min()

    def max(self, ):
        return self.ptr().max()

    def __call__(self, ):
        return deref(self.ptr())()

    def get(self, ):
        return self.ptr().get()

    def getUnfII(self, ):
        return self.ptr().getUnfII()

    def getUnfIE(self, ):
        return self.ptr().getUnfIE()

    def getUnfEE(self, ):
        return self.ptr().getUnfEE()

    def getUnfIE53(self, ):
        return self.ptr().getUnfIE53()

    def getStdExp(self, ):
        return self.ptr().getStdExp()

    def getExp(self, double lambda_):
        return self.ptr().getExp(lambda_)

    def getPsn(self, float lambda_):
        return self.ptr().getPsn(lambda_)

    def getStdNrm(self, ):
        return self.ptr().getStdNrm()

    def getBinom(self, unsigned int t, double p):
        return self.ptr().getBinom(t, p)

    @staticmethod
    cdef _py_RNG from_ptr(RNG *ptr):
        cdef _py_RNG obj = _py_RNG.__new__(_py_RNG)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_RNG from_ref(const RNG &ref):
        return _py_RNG.from_ptr(<RNG*>&ref)


