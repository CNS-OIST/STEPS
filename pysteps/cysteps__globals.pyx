###___license_placeholder___###

from libcpp.string cimport string


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


cdef inline string  to_std_string(str s):
    """
    Builds a C++ STL string from a bytes (Python 2) or unicode (Python 3)
    string.
    """
    cdef string s_new
    if isinstance(s, bytes):
        s_new = s
    else:
        s_new = s.encode()
    return s_new


cdef inline str from_std_string(string s):
    """
    Returns a bytes string (Python 2) or unicode string (Python 3) from a C++
    STL string.
    """
    cdef bytes s_new = s
    return  str(s_new) if bytes is str \
        else s_new.decode()

