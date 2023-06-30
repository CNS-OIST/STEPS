# cython:language_level=3str
###___license_placeholder___###

from libcpp.string cimport string
from libcpp.vector cimport vector

import warnings

cdef enum OPERATOR:
    LESS = 0, LESS_EQUAL, EQUAL, DIFF, GREATER, GREATER_EQUAL


cdef class _py__base:
    cdef void *_ptr

    # Basic comparison is done by comparing the inner obj ptr
    def __richcmp__(_py__base self, _py__base other, operation):
        if operation == OPERATOR.EQUAL:
            return self._ptr==other._ptr

    @property
    def this(self):
        """An identification of the underlying cpp object"""
        cdef char addrstr[20]
        sprintf(addrstr, "%p", self._ptr)
        #string array is converted to a Python2 str (bytes), so that it is ref-counted and is safe to return
        cdef bytes mems = addrstr
        return b"_cPtr_" + mems

    def __hash__(self):
        return hash(<long int>self._ptr)


cdef inline string to_std_string(str s):
    """
    Builds a C++ STL string from a unicode string.
    """
    cdef string s_new
    s_new = s.encode()
    return s_new


cdef inline str from_std_string(string s):
    """
    Returns a unicode string (Python 3) from a C++ STL string.
    """
    cdef bytes s_new = s
    return s_new.decode()

cdef inline vector[string] to_vec_std_strings(str_list):
    """
    Builds a C++ STL vector of STL strings from a list of unicode strings.
    """
    cdef vector[string] str_vec
    str_vec.reserve(len(str_list))
    for s in str_list:
        str_vec.push_back(to_std_string(s))    
    return str_vec

cdef inline list string_set_to_list(const std.set[std.string] &s):
    """
    Return a python list of strings from a C++ STL set of strings.
    """
    return [from_std_string(v) for v in s]

def ShowDeprecationWarning(item, replacement=None, version=None):
    msg = f'{item} is deprecated'
    if version is not None:
        msg += f' and will be removed in STEPS {version}'
    if replacement is not None:
        msg += f', use {replacement} instead.'
    warnings.warn(msg, DeprecationWarning)
