//
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2006 Stefan Wils.
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

////////////////////////////////////////////////////////////////////////////////

// $Id$

////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_MATH_NORMAL3_HPP
#define STEPS_MATH_NORMAL3_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Boost headers.
#include <steps/boost/scoped_ptr.hpp>
#include <steps/boost/shared_ptr.hpp>
#include <steps/boost/weak_ptr.hpp>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/vector3.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

class Normal3;

// Auxiliary declarations.

typedef boost::scoped_ptr<Normal3>              Normal3ScPtr;
typedef boost::shared_ptr<Normal3>              Normal3ShPtr;
typedef boost::weak_ptr<Normal3>                Normal3WkPtr;

////////////////////////////////////////////////////////////////////////////////

class Normal3: 
public Vector3
{

private:

protected:

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Default constructor.
    ///
    Normal3(void)
    {
    }
    
    /// Constructor.
    ///
    Normal3(double x, double y, double z)
    : Vector3(x, y, z)
    {
    }
    
    /// Constructor.
    ///
    Normal3(Vector3 const & v)
    : Vector3(v)
    {
    }
    
    /// Destructor.
    ///
    ~Normal3(void)
    {
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Output operator.
    ///
    friend std::ostream & operator<< (std::ostream & os, Normal3 const & n);
    
};

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN std::ostream & operator<< (std::ostream & os, Normal3 const & n);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_NORMAL3_HPP

// END
