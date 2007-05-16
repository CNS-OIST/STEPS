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

#ifndef STEPS_MATH_GEOMETRY_HPP
#define STEPS_MATH_GEOMETRY_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>

// STEPS headers.
#include <steps/common.h>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

template<typename T>
double area3(T const & v1, T const & v2, T const & v3);

template<typename T>
double area3(T * const vertices[3]);

////////////////////////////////////////////////////////////////////////////////

template<typename T>
double volume4(T const & v1, T const & v2, T const & v3, T const & v4);

template<typename T>
double volume4(T * const vertices[4]);

////////////////////////////////////////////////////////////////////////////////

template<typename T>
void barycenter2(T & c, T const & v1, T const & v2)
{
    c = (v1 + v2) * 0.5;
}

template<typename T>
void barycenter2(T & c, T v[2])
{
    c = (v[0] + v[1]) * 0.5;
}

template<typename T>
void barycenter2(T & c, T * const v[2])
{
    c = (*(v[0]) + *(v[1])) * 0.5;
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
void barycenter3(T & c, T const & v1, T const & v2, T const & v3)
{
    c = (v1 + v2 + v3) * (1.0 / 3.0);
}

template<typename T>
void barycenter3(T & c, T v[3])
{
    c = (v[0] + v[1] + v[2]) * (1.0 / 3.0);
}

template<typename T>
void barycenter3(T & c, T * const v[3])
{
    c = (*(v[0]) + *(v[1]) + *(v[2])) * (1.0 / 3.0); 
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
void barycenter4(T & c, T const & v1, T const & v2, T const & v3, T const & v4)
{
    c = (v1 + v2 + v3 + v4) * 0.25;
}

template<typename T>
void barycenter4(T & c, T v[4])
{
    c = (v[0] + v[1] + v[2] + v[3]) * 0.25;
}

template<typename T>
void barycenter4(T & c, T * const v[4])
{
    c = (*(v[0]) + *(v[1]) + *(v[2]) + *(v[3])) * 0.25;
}

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_GEOMETRY_HPP

// END
