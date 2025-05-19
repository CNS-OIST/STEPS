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
