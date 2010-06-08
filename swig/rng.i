////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006 University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 313 $
 *  Last Changed Date: $Date: 2010-03-25 16:24:21 +0900 (Thu, 25 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

%module rng

%import "cpp/common.h"

%include "python/std_string.i"
%include "error.i"
%{	
#include "../cpp/rng/rng.hpp"
%}


///////////////////////////////////////////////////////////////////////////////	

%exception 
{
	try {
		$action
	} catch (steps::ArgErr & ae) {
		PyErr_SetString(PyExc_NameError, ae.getMsg());
		return NULL;
	} catch (steps::NotImplErr & nie) {
		PyErr_SetString(PyExc_NotImplementedError, nie.getMsg());
		return NULL;
	}
}

///////////////////////////////////////////////////////////////////////////////	

namespace steps
{
	namespace rng
	{
		
		////////////////////////////////////////////////////////////////////////////////
		
        %feature("autodoc", 
"
Proxy of C++ random number generator.

");	
		class RNG
		{
			
		public:
			
			RNG(unsigned int bufsize);
			virtual ~RNG(void);
			
            %feature("autodoc", 
"
Initialize the random number generator with given seed value.


Syntax::
    
    initialize(seed)
    
Arguments:
    uint seed

Return:
    None
");
			void initialize(unsigned long const & seed);
			
			unsigned int get(void);
			
			double getUnfII(void);
			double getUnfIE(void);
			double getUnfEE(void);
			double getUnfIE53(void);
			
			float getStdExp(void);
			virtual double getExp(double lambda);
			long getPsn(float lambda);
			float getStdNrm(void);
			
		protected:
			// Mark this class as abstract.
			virtual void concreteInitialize(ulong seed) = 0;
			virtual void concreteFillBuffer(void) = 0;
			
		};
		
		////////////////////////////////////////////////////////////////////////////////
        %feature("autodoc", 
"
Equivalent to: create('mt19937', buffer_size)

Syntax::

    create_mt19937(buffer_size)
    
Arguments:
    uint buffer_size

Return:
    steps.rng.RNG
");		
		RNG * create_mt19937(unsigned int bufsize);
		
        %feature("autodoc", 
"
Creates and returns a reference to a steps.rng.RNG random number generator object, 
which is specified by type and pre-allocates a buffer list with size of buffer_size.

Syntax::
    
    create(type, buffer_size)

Arguments:
    * string type
    * uint buffer_size

Return:
    steps.rng.RNG
");	
		RNG * create(std::string rng_name, unsigned int bufsize);
		
		////////////////////////////////////////////////////////////////////////////////
		
	} // end namespace rng
} // end namespace steps

// END

