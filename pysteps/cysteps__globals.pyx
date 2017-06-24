###___license_placeholder___###

ctypedef unsigned int uint
ctypedef unsigned short ushort
ctypedef unsigned long ulong
ctypedef unsigned char uchar
ctypedef unsigned int size_t


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
        return "_cPtr_" + mems

