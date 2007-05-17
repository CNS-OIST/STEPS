////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <cmath>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/geometry.hpp>
#include <steps/math/matrix4.hpp>
#include <steps/math/point3.hpp>
#include <steps/math/vector3.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
USING(smath, Matrix4);
USING(smath, Point3);
USING(smath, Vector3);

////////////////////////////////////////////////////////////////////////////////

// Generic implementation.

template<typename T>
double smath::area3(T const & v1, T const & v2, T const & v3)
{
    return 0.5 * ((v2 - v1) ^ (v3 - v1)).getLength();
}

// Explicit instantiation.

template double smath::area3<Point3>
(
    Point3 const & v1, 
    Point3 const & v2, 
    Point3 const & v3
);

template double smath::area3<Vector3>
(
    Vector3 const & v1, 
    Vector3 const & v2, 
    Vector3 const & v3
);

////////////////////////////////////////////////////////////////////////////////

// Generic implementation.

template<typename T>
double smath::area3(T * const v[3])
{
    return 0.5 * ((*(v[1]) - *(v[0])) ^ (*(v[2]) - *(v[0]))).getLength();
}

// Explicit instantiation.

template double smath::area3<Point3>(Point3 * const v[3]);

template double smath::area3<Vector3>(Vector3 * const v[3]);

////////////////////////////////////////////////////////////////////////////////

// Generic implementation.

template<typename T>
double smath::volume4(T const & v1, T const & v2, T const & v3, T const & v4)
{
    Matrix4 m(v1.getX(), v1.getY(), v1.getZ(), 1.0,
              v2.getX(), v2.getY(), v2.getZ(), 1.0,
              v3.getX(), v3.getY(), v3.getZ(), 1.0,
              v4.getX(), v4.getY(), v4.getZ(), 1.0);
    return fabs(m.determinant() / 6.0);
}

// Explicit instantiation.

template double smath::volume4<Point3>
(
    Point3 const & v1, 
    Point3 const & v2, 
    Point3 const & v3, 
    Point3 const & v4
);

template double smath::volume4<Vector3>
(
    Vector3 const & v1, 
    Vector3 const & v2, 
    Vector3 const & v3, 
    Vector3 const & v4
);

////////////////////////////////////////////////////////////////////////////////

// Generic implementation.

template<typename T>
double smath::volume4(T * const v[4])
{
    Matrix4 m(v[0]->getX(), v[0]->getY(), v[0]->getZ(), 1.0,
              v[1]->getX(), v[1]->getY(), v[1]->getZ(), 1.0,
              v[2]->getX(), v[2]->getY(), v[2]->getZ(), 1.0,
              v[3]->getX(), v[3]->getY(), v[3]->getZ(), 1.0);
    m.transpose();
    return fabs(m.determinant() / 6.0);
}

// Explicit instantiation.

template double smath::volume4<Point3>(Point3 * const v[4]);

template double smath::volume4<Vector3>(Vector3 * const v[4]);

////////////////////////////////////////////////////////////////////////////////

// END
