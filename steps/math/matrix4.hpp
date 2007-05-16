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

#ifndef STEPS_MATH_MATRIX4_HPP
#define STEPS_MATH_MATRIX4_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <cmath>
#include <iostream>

// Boost headers.
#include <steps/boost/scoped_ptr.hpp>
#include <steps/boost/shared_ptr.hpp>
#include <steps/boost/weak_ptr.hpp>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/aabbox3.hpp>
#include <steps/math/normal3.hpp>
#include <steps/math/point3.hpp>
#include <steps/math/tools.hpp>
#include <steps/math/vector3.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

class Matrix4;

typedef boost::scoped_ptr<Matrix4>              Matrix4ScPtr;
typedef boost::shared_ptr<Matrix4>              Matrix4ShPtr;
typedef boost::weak_ptr<Matrix4>                Matrix4WkPtr;

////////////////////////////////////////////////////////////////////////////////

class Matrix4
{

private:

protected:

    double                                      rC[4][4];

public:

    Matrix4(void)
    {
    }
    
    Matrix4(double * a)
    {
        rC[0][0] = a[ 0]; rC[0][1] = a[ 1]; rC[0][2] = a[ 2]; rC[0][3] = a[ 3];
        rC[1][0] = a[ 4]; rC[1][1] = a[ 5]; rC[1][2] = a[ 6]; rC[1][3] = a[ 7];
        rC[2][0] = a[ 8]; rC[2][1] = a[ 9]; rC[2][2] = a[10]; rC[2][3] = a[11];
        rC[3][0] = a[12]; rC[3][1] = a[13]; rC[3][2] = a[14]; rC[3][3] = a[15];
    }
    
    Matrix4(double a11, double a12, double a13, double a14,
            double a21, double a22, double a23, double a24,
            double a31, double a32, double a33, double a34,
            double a41, double a42, double a43, double a44)
    {
        rC[0][0] = a11; rC[0][1] = a12; rC[0][2] = a13; rC[0][3] = a14;
        rC[1][0] = a21; rC[1][1] = a22; rC[1][2] = a23; rC[1][3] = a24;
        rC[2][0] = a31; rC[2][1] = a32; rC[2][2] = a33; rC[2][3] = a34;
        rC[3][0] = a41; rC[3][1] = a42; rC[3][2] = a43; rC[3][3] = a44;
    }
    
    ~Matrix4(void)
    {
    }

    static Matrix4 Identity;
    
    double get(int r, int c) const
    {
        assert(r >= 0);
        assert(r < 4);
        assert(c >= 0);
        assert(c < 4);
        return rC[r][c];
    }
    
    double get11(void) const
    {
        return rC[0][0];
    }
    
    double get12(void) const
    {
        return rC[0][1];
    }
    
    double get13(void) const
    {
        return rC[0][2];
    }
    
    double get14(void) const
    {
        return rC[0][3];
    }
    
    double get21(void) const
    {
        return rC[1][0];
    }
    
    double get22(void) const
    {
        return rC[1][1];
    }
    
    double get23(void) const
    {
        return rC[1][2];
    }
    
    double get24(void) const
    {
        return rC[1][3];
    }
    
    double get31(void) const
    {
        return rC[2][0];
    }
    
    double get32(void) const
    {
        return rC[2][1];
    }
    
    double get33(void) const
    {
        return rC[2][2];
    }
    
    double get34(void) const
    {
        return rC[2][3];
    }
    
    double get41(void) const
    {
        return rC[3][0];
    }
    
    double get42(void) const
    {
        return rC[3][1];
    }
    
    double get43(void) const
    {
        return rC[3][2];
    }
    
    double get44(void) const
    {
        return rC[3][3];
    }

    void set(int r, int c, double v)
    {
        assert(r >= 0);
        assert(r < 4);
        assert(c >= 0);
        assert(c < 4);
        rC[r][c] = v;
    }
    
    void set(double a11, double a12, double a13, double a14,
             double a21, double a22, double a23, double a24,
             double a31, double a32, double a33, double a34,
             double a41, double a42, double a43, double a44)
    {
        rC[0][0] = a11; rC[0][1] = a12; rC[0][2] = a13; rC[0][3] = a14;
        rC[1][0] = a21; rC[1][1] = a22; rC[1][2] = a23; rC[1][3] = a24;
        rC[2][0] = a31; rC[2][1] = a32; rC[2][2] = a33; rC[2][3] = a34;
        rC[3][0] = a41; rC[3][1] = a42; rC[3][2] = a43; rC[3][3] = a44;
    }
    
    void set11(double v)
    {
        rC[0][0] = v;
    }
    
    void set12(double v)
    {
        rC[0][1] = v;
    }
    
    void set13(double v)
    {
        rC[0][2] = v;
    }
    
    void set14(double v)
    {
        rC[0][3] = v;
    }
    
    void set21(double v)
    {
        rC[1][0] = v;
    }
    
    void set22(double v)
    {
        rC[1][1] = v;
    }
    
    void set23(double v)
    {
        rC[1][2] = v;
    }
    
    void set24(double v)
    {
        rC[1][3] = v;
    }
    
    void set31(double v)
    {
        rC[2][0] = v;
    }
    
    void set32(double v)
    {
        rC[2][1] = v;
    }
    
    void set33(double v)
    {
        rC[2][2] = v;
    }
    
    void set34(double v)
    {
        rC[2][3] = v;
    }
    
    void set41(double v)
    {
        rC[3][0] = v;
    }
    
    void set42(double v)
    {
        rC[3][1] = v;
    }
    
    void set43(double v)
    {
        rC[3][2] = v;
    }
    
    void set44(double v)
    {
        rC[3][3] = v;
    }
    
    bool isAlmostEqual(const Matrix4& m, double tolerance) const
    {
        assert(tolerance >= 0.0f);
        if (fabsf(rC[0][0] - m.rC[0][0]) > tolerance) return false;
        if (fabsf(rC[0][1] - m.rC[0][1]) > tolerance) return false;
        if (fabsf(rC[0][2] - m.rC[0][2]) > tolerance) return false;
        if (fabsf(rC[0][3] - m.rC[0][3]) > tolerance) return false;
        if (fabsf(rC[1][0] - m.rC[1][0]) > tolerance) return false;
        if (fabsf(rC[1][1] - m.rC[1][1]) > tolerance) return false;
        if (fabsf(rC[1][2] - m.rC[1][2]) > tolerance) return false;
        if (fabsf(rC[1][3] - m.rC[1][3]) > tolerance) return false;
        if (fabsf(rC[2][0] - m.rC[2][0]) > tolerance) return false;
        if (fabsf(rC[2][1] - m.rC[2][1]) > tolerance) return false;
        if (fabsf(rC[2][2] - m.rC[2][2]) > tolerance) return false;
        if (fabsf(rC[2][3] - m.rC[2][3]) > tolerance) return false;
        if (fabsf(rC[3][0] - m.rC[3][0]) > tolerance) return false;
        if (fabsf(rC[3][1] - m.rC[3][1]) > tolerance) return false;
        if (fabsf(rC[3][2] - m.rC[3][2]) > tolerance) return false;
        if (fabsf(rC[3][3] - m.rC[3][3]) > tolerance) return false;
        return true;
    }
    
    void setIdentity(void);
    
    void setRotateX(double a);
    
    void setRotateY(double a);
    
    void setRotateZ(double a);
    
    void setScale(double s);
    
    void setScale(double x, double y, double z);
    
    void setScale(Vector3 const & s);
    
    void setTranslate(double x, double y, double z);
    
    void setTranslate(Vector3 const & t);
    
    double determinant(void) const;
    
    /// \todo       {The current algorithm is Gauss-Jordan with partial
    ///              pivoting. Its stability is sub-optimal, maybe use
    ///              another algorithm?}
    void invert(void);
    
    void transpose(void);
    
    friend bool operator== (Matrix4 const & m1, Matrix4 const & m2)
    {
        if (steps::math::isLargerEps(m1.rC[0][0] - m2.rC[0][0])) return false;
        if (steps::math::isLargerEps(m1.rC[0][1] - m2.rC[0][1])) return false;
        if (steps::math::isLargerEps(m1.rC[0][2] - m2.rC[0][2])) return false;
        if (steps::math::isLargerEps(m1.rC[0][3] - m2.rC[0][3])) return false;
        if (steps::math::isLargerEps(m1.rC[1][0] - m2.rC[1][0])) return false;
        if (steps::math::isLargerEps(m1.rC[1][1] - m2.rC[1][1])) return false;
        if (steps::math::isLargerEps(m1.rC[1][2] - m2.rC[1][2])) return false;
        if (steps::math::isLargerEps(m1.rC[1][3] - m2.rC[1][3])) return false;
        if (steps::math::isLargerEps(m1.rC[2][0] - m2.rC[2][0])) return false;
        if (steps::math::isLargerEps(m1.rC[2][1] - m2.rC[2][1])) return false;
        if (steps::math::isLargerEps(m1.rC[2][2] - m2.rC[2][2])) return false;
        if (steps::math::isLargerEps(m1.rC[2][3] - m2.rC[2][3])) return false;
        if (steps::math::isLargerEps(m1.rC[3][0] - m2.rC[3][0])) return false;
        if (steps::math::isLargerEps(m1.rC[3][1] - m2.rC[3][1])) return false;
        if (steps::math::isLargerEps(m1.rC[3][2] - m2.rC[3][2])) return false;
        if (steps::math::isLargerEps(m1.rC[3][3] - m2.rC[3][3])) return false;
        return true;
    }
    
    friend bool operator!= (Matrix4 const & m1, Matrix4 const & m2)
    {
        if (steps::math::isLargerEps(m1.rC[0][0] - m2.rC[0][0])) return true;
        if (steps::math::isLargerEps(m1.rC[0][1] - m2.rC[0][1])) return true;
        if (steps::math::isLargerEps(m1.rC[0][2] - m2.rC[0][2])) return true;
        if (steps::math::isLargerEps(m1.rC[0][3] - m2.rC[0][3])) return true;
        if (steps::math::isLargerEps(m1.rC[1][0] - m2.rC[1][0])) return true;
        if (steps::math::isLargerEps(m1.rC[1][1] - m2.rC[1][1])) return true;
        if (steps::math::isLargerEps(m1.rC[1][2] - m2.rC[1][2])) return true;
        if (steps::math::isLargerEps(m1.rC[1][3] - m2.rC[1][3])) return true;
        if (steps::math::isLargerEps(m1.rC[2][0] - m2.rC[2][0])) return true;
        if (steps::math::isLargerEps(m1.rC[2][1] - m2.rC[2][1])) return true;
        if (steps::math::isLargerEps(m1.rC[2][2] - m2.rC[2][2])) return true;
        if (steps::math::isLargerEps(m1.rC[2][3] - m2.rC[2][3])) return true;
        if (steps::math::isLargerEps(m1.rC[3][0] - m2.rC[3][0])) return true;
        if (steps::math::isLargerEps(m1.rC[3][1] - m2.rC[3][1])) return true;
        if (steps::math::isLargerEps(m1.rC[3][2] - m2.rC[3][2])) return true;
        if (steps::math::isLargerEps(m1.rC[3][3] - m2.rC[3][3])) return true;
        return false;
    }

    friend AABBox3 operator* (Matrix4 const & m, AABBox3 const & bb);

    friend Matrix4 operator* (Matrix4 const & m1, Matrix4 const & m2);
    
    /// \todo       {Get rid of the transposition step, because:
    ///              (M * v) == (v^T * M^T). Also, maybe a shortcut for 
    ///              the inversion should be provided, as some extra 
    ///              xformNormal() function that assumes that the matrix 
    ///              is already inverted.}
    friend Normal3 operator* (Matrix4 const & m, Normal3 const & n);
    
    friend Point3 operator* (Matrix4 const & m, Point3 const & p);
    
    friend Vector3 operator* (Matrix4 const & m, Vector3 const & v);

    ////////////////////////////////////////////////////////////////////////

    friend std::istream & operator>> (std::istream & is, Matrix4 & m);

    friend std::ostream & operator<< (std::ostream & os, Matrix4 const & m);
    
};

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN AABBox3 operator* (Matrix4 const & m, AABBox3 const & bb);

STEPS_EXTERN Matrix4 operator* (Matrix4 const & m1, Matrix4 const & m2);

STEPS_EXTERN Normal3 operator* (Matrix4 const & m, Normal3 const & n);

STEPS_EXTERN Point3 operator* (Matrix4 const & m, Point3 const & p);

STEPS_EXTERN Vector3 operator* (Matrix4 const & m, Vector3 const & v);

STEPS_EXTERN std::istream & operator>> (std::istream & is, Matrix4 & m);

STEPS_EXTERN std::ostream & operator<< (std::ostream & os, Matrix4 const & m);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_MATRIX4_HPP

// END
