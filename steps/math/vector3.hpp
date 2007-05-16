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

#ifndef STEPS_MATH_VECTOR3_HPP
#define STEPS_MATH_VECTOR3_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cmath>
#include <iostream>

// Boost headers.
#include <steps/boost/scoped_ptr.hpp>
#include <steps/boost/shared_ptr.hpp>
#include <steps/boost/weak_ptr.hpp>

// STEPS headers.
#include <steps/math/tools.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

class Vector3;

// Auxiliary declarations.

typedef boost::scoped_ptr<Vector3>              Vector3ScPtr;
typedef boost::shared_ptr<Vector3>              Vector3ShPtr;
typedef boost::weak_ptr<Vector3>                Vector3WkPtr;

////////////////////////////////////////////////////////////////////////////////

/// \todo                   {Add an input stream operator ">>" to read vectors,
///                         normals and points from a general character-based
///                         input stream.}

class Vector3
{

private:

protected:

    double                                      rX;
    
    double                                      rY;
    
    double                                      rZ;

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Default constructor.
    ///
    Vector3(void)
    {
    }
    
    /// Constructor.
    ///
    Vector3(double x, double y, double z)
    : rX(x)
    , rY(y)
    , rZ(z)
    {
    }

    /// Destructor.
    ///
    ~Vector3(void)
    {
    }
    
    ////////////////////////////////////////////////////////////////////////
    // VECTOR COMPONENT ACCESS
    ////////////////////////////////////////////////////////////////////////
    
    double getX(void) const
    {
        return rX;
    }
    
    double getY(void) const
    {
        return rY;
    }
    
    double getZ(void) const
    {
        return rZ;
    }

    void setX(double x)
    {
        rX = x;
    }
    
    void setY(double y)
    {
        rY = y;
    }
    
    void setZ(double z)
    {
        rZ = z;
    }

    void set(double x, double y, double z)
    {
        rX = x;
        rY = y;
        rZ = z;
    }
    
    bool isAlmostEqual(Vector3 const & v, double tolerance) const
    {
        assert(tolerance >= 0.0);
        if (fabsf(rX - v.rX) > tolerance) return false;
        if (fabsf(rY - v.rY) > tolerance) return false;
        if (fabsf(rZ - v.rZ) > tolerance) return false;
        return true;
    }
    
    double getLength(void) const
    {
        return sqrt((rX * rX) + (rY * rY) + (rZ * rZ));
    }
    
    double getSquaredLength(void) const
    {
        return (rX * rX) + (rY * rY) + (rZ * rZ);
    }
    
    void setNormalize(void);

    Vector3 & operator= (float v[3])
    {
        rX = v[0];
        rY = v[1];
        rZ = v[2];
        return (*this);
    }

    Vector3 & operator= (double v[3])
    {
        rX = v[0];
        rY = v[1];
        rZ = v[2];
        return (*this);
    }

    Vector3 & operator+= (Vector3 const & v)
    {
        rX += v.rX;
        rY += v.rY;
        rZ += v.rZ;
        return (*this);
    }
    
    Vector3 & operator+= (double r)
    {
        rX += r;
        rY += r;
        rZ += r;
        return (*this);
    }

    Vector3 & operator-= (Vector3 const & v)
    {
        rX -= v.rX;
        rY -= v.rY;
        rZ -= v.rZ;
        return (*this);
    }
    
    Vector3 & operator-= (double r)
    {
        rX -= r;
        rY -= r;
        rZ -= r;
        return (*this);
    }
    
    Vector3 & operator*= (Vector3 const & v)
    {
        rX *= v.rX;
        rY *= v.rY;
        rZ *= v.rZ;
        return (*this);
    }
    
    Vector3 & operator*= (double r)
    {
        rX *= r;
        rY *= r;
        rZ *= r;
        return (*this);
    }

    Vector3 & operator/= (Vector3 const & v)
    {
        rX /= v.rX;
        rY /= v.rY;
        rZ /= v.rZ;
        return (*this);
    }
    
    Vector3 & operator/= (double r)
    {
        rX /= r;
        rY /= r;
        rZ /= r;
        return (*this);
    }
    
    friend bool operator== (Vector3 const & v1, Vector3 const & v2)
    {
        if (!steps::math::isSmallerEps(v1.rX - v2.rX)) return false;
        if (!steps::math::isSmallerEps(v1.rY - v2.rY)) return false;
        if (!steps::math::isSmallerEps(v1.rZ - v2.rZ)) return false;
        return true;
    }
    
    friend bool operator!= (Vector3 const & v1, Vector3 const & v2)
    {
        return (!steps::math::isSmallerEps(v1.rX - v2.rX))
            || (!steps::math::isSmallerEps(v1.rY - v2.rY))
            || (!steps::math::isSmallerEps(v1.rZ - v2.rZ));
    }
    
    friend Vector3 operator+ (Vector3 const & v1, Vector3 const & v2)
    {
        Vector3 r(v1);
        r += v2;
        return r;
    }
        
    friend Vector3 operator- (Vector3 const & v)
    {
        return Vector3(-v.rX, -v.rY, -v.rZ);
    }
    
    friend Vector3 operator- (Vector3 const & v1, Vector3 const & v2)
    {
        Vector3 r(v1);
        r -= v2;
        return r;
    }
    
    friend Vector3 operator* (Vector3 const & v, double s)
    {
        Vector3 r(v);
        r *= s;
        return r;
    }
    
    friend Vector3 operator* (double s, Vector3 const & v)
    {
        Vector3 r(v);
        r *= s;
        return r;
    }
    
    friend Vector3 operator/ (Vector3 const & v1, Vector3 const & v2)
    {
        Vector3 r(v1);
        r /= v2;
        return r;
    }
    
    friend double operator* (Vector3 const & v1, Vector3 const & v2)
    {
        return (v1.rX * v2.rX) + (v1.rY * v2.rY) + (v1.rZ * v2.rZ);
    }

    friend Vector3 operator^ (Vector3 const & v1, Vector3 const & v2);
    
    ////////////////////////////////////////////////////////////////////////
    
    friend std::istream & operator>> (std::istream & is, Vector3 & vec);
    
    friend std::ostream & operator<< (std::ostream & os, Vector3 const & vec);
    
};

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN Vector3 operator^ (Vector3 const & v1, Vector3 const & v2);

STEPS_EXTERN std::istream & operator>> (std::istream & is, Vector3 & vec);

STEPS_EXTERN std::ostream & operator<< (std::ostream & os, Vector3 const & vec);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_VECTOR3_HPP

// END
