###___license_placeholder___###

cdef extern from "util/common.hpp":
    ctypedef unsigned int uint
    ctypedef unsigned short ushort
    ctypedef unsigned long ulong
    ctypedef unsigned char uchar
    ctypedef unsigned int size_t

cdef extern from "util/vocabulary.hpp" namespace "steps":
    ctypedef char index_t;

cdef extern from "util/collections.hpp" namespace "steps::util" nogil:
   cdef cppclass flat_set[T]:
       cppclass iterator:
           T& operator*()
           iterator operator++()
           iterator operator--()
           iterator operator+(size_type)
           iterator operator-(size_type)
           bint operator==(iterator)
           bint operator!=(iterator)
           bint operator<(iterator)
           bint operator>(iterator)
           bint operator<=(iterator)
           bint operator>=(iterator)

       cppclass const_iterator(iterator):
           pass

       flat_set() except +
       iterator begin()
       const_iterator const_begin "begin"()
       const_iterator const_end "end"()
       iterator end()
