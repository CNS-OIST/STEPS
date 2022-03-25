###___license_placeholder___###

cdef extern from "util/common.h":
    ctypedef unsigned int uint
    ctypedef unsigned short ushort
    ctypedef unsigned long ulong
    ctypedef unsigned char uchar
    ctypedef unsigned int size_t

cdef extern from "util/vocabulary.hpp" namespace "steps":
    ctypedef char index_t;
