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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/vector3.hpp>
#include <steps/math/tools.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
USING(smath, Vector3);

////////////////////////////////////////////////////////////////////////////////

void Vector3::setNormalize(void)
{
    // Compute the length.
    double l = sqrt((rX * rX) + (rY * rY) + (rZ * rZ));
    // Check whether it is positive.
    assert(l > 0.0);
    // Normalize.
    rX /= l;
    rY /= l;
    rZ /= l;
}

////////////////////////////////////////////////////////////////////////////////

Vector3 smath::operator^ (Vector3 const & v1, Vector3 const & v2)
{
    return (Vector3((v1.rY * v2.rZ) - (v1.rZ * v2.rY),
                    (v1.rZ * v2.rX) - (v1.rX * v2.rZ),
                    (v1.rX * v2.rY) - (v1.rY * v2.rX)));
}

////////////////////////////////////////////////////////////////////////////////

std::istream & smath::operator>> (std::istream & is, Vector3 & vec)
{
    // Buffer to store the components; vec should remain untouched if the
    // read operation fails.
    double x, y, z;
    // Initialize c; if >> fails, c should not be set!
    char c = 0;
    
    // Read (skipping whitespace).
    is >> c;
    if (c != '(') goto error;
    is >> x >> c;
    if (c != ',') goto error;
    is >> y >> c;
    if (c != ',') goto error;
    is >> z >> c;
    if (c != ')') goto error;
    
    // Copy the vector.
    vec = Vector3(x, y, z);
    
final:
    return is;
    
error:
    is.clear(std::ios_base::badbit);
    goto final;
}

////////////////////////////////////////////////////////////////////////////////

std::ostream & smath::operator<< (std::ostream & os, Vector3 const & vec)
{
    os << "(" << vec.getX() << ", "
              << vec.getY() << ", "
              << vec.getZ() << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////

// END
