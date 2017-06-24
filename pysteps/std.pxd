from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string
from libcpp.queue cimport queue
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport pair

cdef extern from "exception" namespace "std":
    cdef cppclass exception:
        pass

cdef extern from "stdexcept" namespace "std":
    cdef cppclass runtime_error:
        pass
