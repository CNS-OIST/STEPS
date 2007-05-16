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
#include <cmath>
#include <cstring>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/aabbox3.hpp>
#include <steps/math/matrix4.hpp>
#include <steps/math/normal3.hpp>
#include <steps/math/point3.hpp>
#include <steps/math/tools.hpp>
#include <steps/math/vector3.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
USING(smath, AABBox3);
USING(smath, Matrix4);
USING(smath, Normal3);
USING(smath, Point3);
USING(smath, Vector3);

////////////////////////////////////////////////////////////////////////////////

Matrix4 Matrix4::Identity
(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
);

////////////////////////////////////////////////////////////////////////////////

void Matrix4::setIdentity(void)
{
    rC[0][0] = 1.0f; rC[0][1] = 0.0f; rC[0][2] = 0.0f; rC[0][3] = 0.0f;
    rC[1][0] = 0.0f; rC[1][1] = 1.0f; rC[1][2] = 0.0f; rC[1][3] = 0.0f;
    rC[2][0] = 0.0f; rC[2][1] = 0.0f; rC[2][2] = 1.0f; rC[2][3] = 0.0f;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setRotateX(double a)
{
    double sa = sin(a);
    double ca = cos(a);
    rC[0][0] = 1.0f; rC[0][1] = 0.0f; rC[0][2] = 0.0f; rC[0][3] = 0.0f;
    rC[1][0] = 0.0f; rC[1][1] = ca;   rC[1][2] = -sa;  rC[1][3] = 0.0f;
    rC[2][0] = 0.0f; rC[2][1] = sa;   rC[2][2] = ca;   rC[2][3] = 0.0f;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setRotateY(double a)
{
    double sa = sin(a);
    double ca = cos(a);
    rC[0][0] = ca;   rC[0][1] = 0.0f; rC[0][2] = sa;   rC[0][3] = 0.0f;
    rC[1][0] = 0.0f; rC[1][1] = 1.0f; rC[1][2] = 0.0f; rC[1][3] = 0.0f;
    rC[2][0] = -sa;  rC[2][1] = 0.0f; rC[2][2] = ca;   rC[2][3] = 0.0f;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setRotateZ(double a)
{
    double sa = sin(a);
    double ca = cos(a);
    rC[0][0] = ca;   rC[0][1] = -sa;  rC[0][2] = 0.0f; rC[0][3] = 0.0f;
    rC[1][0] = sa;   rC[1][1] = ca;   rC[1][2] = 0.0f; rC[1][3] = 0.0f;
    rC[2][0] = 0.0f; rC[2][1] = 0.0f; rC[2][2] = 1.0f; rC[2][3] = 0.0f;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setScale(double s)
{
    rC[0][0] = s;    rC[0][1] = 0.0f; rC[0][2] = 0.0f; rC[0][3] = 0.0f;
    rC[1][0] = 0.0f; rC[1][1] = s;    rC[1][2] = 0.0f; rC[1][3] = 0.0f;
    rC[2][0] = 0.0f; rC[2][1] = 0.0f; rC[2][2] = s;    rC[2][3] = 0.0f;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////    
    
void Matrix4::setScale(double x, double y, double z)
{
    rC[0][0] = x;    rC[0][1] = 0.0f; rC[0][2] = 0.0f; rC[0][3] = 0.0f;
    rC[1][0] = 0.0f; rC[1][1] = y;    rC[1][2] = 0.0f; rC[1][3] = 0.0f;
    rC[2][0] = 0.0f; rC[2][1] = 0.0f; rC[2][2] = z;    rC[2][3] = 0.0f;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setScale(Vector3 const & s)
{
    rC[0][0] = s.getX(); rC[0][1] = 0.0f;
    rC[0][2] = 0.0f;     rC[0][3] = 0.0f;
    rC[1][0] = 0.0f;     rC[1][1] = s.getY();
    rC[1][2] = 0.0f;     rC[1][3] = 0.0f;
    rC[2][0] = 0.0f;     rC[2][1] = 0.0f;
    rC[2][2] = s.getZ(); rC[2][3] = 0.0f;
    rC[3][0] = 0.0f;     rC[3][1] = 0.0f;
    rC[3][2] = 0.0f;     rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setTranslate(double x, double y, double z)
{
    rC[0][0] = 1.0f; rC[0][1] = 0.0f; rC[0][2] = 0.0f; rC[0][3] = x;
    rC[1][0] = 0.0f; rC[1][1] = 1.0f; rC[1][2] = 0.0f; rC[1][3] = y;
    rC[2][0] = 0.0f; rC[2][1] = 0.0f; rC[2][2] = 1.0f; rC[2][3] = z;
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::setTranslate(Vector3 const & t)
{
    rC[0][0] = 1.0f; rC[0][1] = 0.0f; rC[0][2] = 0.0f; rC[0][3] = t.getX();
    rC[1][0] = 0.0f; rC[1][1] = 1.0f; rC[1][2] = 0.0f; rC[1][3] = t.getY();
    rC[2][0] = 0.0f; rC[2][1] = 0.0f; rC[2][2] = 1.0f; rC[2][3] = t.getZ();
    rC[3][0] = 0.0f; rC[3][1] = 0.0f; rC[3][2] = 0.0f; rC[3][3] = 1.0f;
}
    
////////////////////////////////////////////////////////////////////////////////

inline double det3x3(double a1, double a2, double a3, 
                    double b1, double b2, double b3,
                    double c1, double c2, double c3)
{
    return a1 * ((b2 * c3) - (c2 * b3))
         - b1 * ((a2 * c3) - (c2 * a3))
         + c1 * ((a2 * b3) - (b2 * a3));
}

double Matrix4::determinant(void) const
{
    double a30 = rC[3][0];
    double a31 = rC[3][1];
    double a32 = rC[3][2];
    
    if (smath::isSmallerEps(a30) &&
        smath::isSmallerEps(a31) && 
        smath::isSmallerEps(a32))
        return -rC[3][0] * det3x3(rC[0][1], rC[0][2], rC[0][3],
                                  rC[1][1], rC[1][2], rC[1][3],
                                  rC[2][1], rC[2][2], rC[2][3]);
    
    return rC[0][0] * (det3x3(rC[1][1], rC[1][2], rC[1][3],
                              rC[2][1], rC[2][2], rC[2][3],
                              a31, a32, rC[3][3]))
         - rC[1][0] * (det3x3(rC[0][1], rC[0][2], rC[0][3],
                              rC[2][1], rC[2][2], rC[2][3],
                              a31, a32, rC[3][3]))
         + rC[2][0] * (det3x3(rC[0][1], rC[0][2], rC[0][3],
                              rC[1][1], rC[1][2], rC[1][3],
                              a31, a32, rC[3][3]))
         - a30 *      (det3x3(rC[0][1], rC[0][2], rC[0][3],
                              rC[1][1], rC[1][2], rC[1][3],
                              rC[2][1], rC[2][2], rC[2][3]));
}

////////////////////////////////////////////////////////////////////////////////

void Matrix4::invert(void)
{
    // Create a temporary matrix that will store the inverse.
    // TODO: maybe we can copy some pre-made identity matrix for speed?
    double iv[4][4];
    iv[0][0] = 1.0f; iv[0][1] = 0.0f; iv[0][2] = 0.0f; iv[0][3] = 0.0f;
    iv[1][0] = 0.0f; iv[1][1] = 1.0f; iv[1][2] = 0.0f; iv[1][3] = 0.0f;
    iv[2][0] = 0.0f; iv[2][1] = 0.0f; iv[2][2] = 1.0f; iv[2][3] = 0.0f;
    iv[3][0] = 0.0f; iv[3][1] = 0.0f; iv[3][2] = 0.0f; iv[3][3] = 1.0f;
    
    // Loop over all 4 columns of the matrix.
    for (int i = 0; i < 4; i++) {
    
        // PARTIAL PIVOTING: select a pivot for the i-th column, store
        // it in k. If no pivot can be found, we have a singular matrix.
        // Else, we swap the pivot row with the current (i-th) row.
        
        // pvt stores the pivot's row index.
        int pvt = i;
        // pvt_v stores the absolute value of the pivot.
        double pvt_v = fabs(rC[pvt][pvt]);
        // pvt_t does temporary storage.
        double pvt_t;
        
        switch (i) {
        case 0:
            pvt_t = fabs(rC[1][i]);
            if (pvt_t > pvt_v) {
                pvt_v = pvt_t;
                pvt = 1;
            }
        case 1:
            pvt_t = fabs(rC[2][i]);
            if (pvt_t > pvt_v) {
                pvt_v = pvt_t;
                pvt = 2;
            }
        case 2:
            pvt_t = fabs(rC[3][i]);
            if (pvt_t > pvt_v) {
                pvt_v = pvt_t;
                pvt = 3;
            }
        }
        
        // Do we have a singular matrix?
        assert(smath::isSmallerEps(pvt_v) == false);
        
        // Else, swap the pvt-th row with the i-th row.
        if (i != pvt) {
            double buf;
            switch (i) {
            case 0:
                buf = rC[i][0]; rC[i][0] = rC[pvt][0]; rC[pvt][0] = buf;
            case 1:
                buf = rC[i][1]; rC[i][1] = rC[pvt][1]; rC[pvt][1] = buf;
            case 2:
                buf = rC[i][2]; rC[i][2] = rC[pvt][2]; rC[pvt][2] = buf;
            case 3:
                buf = rC[i][3]; rC[i][3] = rC[pvt][3]; rC[pvt][3] = buf;
            }
            buf = iv[i][0]; iv[i][0] = iv[pvt][0]; iv[pvt][0] = buf;
            buf = iv[i][1]; iv[i][1] = iv[pvt][1]; iv[pvt][1] = buf;
            buf = iv[i][2]; iv[i][2] = iv[pvt][2]; iv[pvt][2] = buf;
            buf = iv[i][3]; iv[i][3] = iv[pvt][3]; iv[pvt][3] = buf;
        }
        
        // ELIMINATION: the pivot row should be scaled, all other rows must
        // have the element that is in the pivot column zero-ed out, by making
        // a linear combination of the row with the pivot row.

        // Fetch the pivot value.
        pvt_v = rC[i][i];
        
        // Scale the row with the pivot value.
        switch (i) {
        case 0:
            rC[i][1] /= pvt_v;
        case 1:
            rC[i][2] /= pvt_v;
        case 2:
            rC[i][3] /= pvt_v;
        }
        iv[i][0] /= pvt_v;
        iv[i][1] /= pvt_v;
        iv[i][2] /= pvt_v;
        iv[i][3] /= pvt_v;
        
        // Zero out the entire i-th column: make linear combinations of rows.
        for (int j = 0; j < 4; j++) {
            if (j != i) {
                double x = rC[j][i];
                switch (i) {
                case 0:
                    rC[j][1] -= rC[i][1] * x;
                case 1:
                    rC[j][2] -= rC[i][2] * x;
                case 2:
                    rC[j][3] -= rC[i][3] * x;
                }
                iv[j][0] -= iv[i][0] * x;
                iv[j][1] -= iv[i][1] * x;
                iv[j][2] -= iv[i][2] * x;
                iv[j][3] -= iv[i][3] * x;
            }
        }
    }
    
    // Copy the temporary matrix into the object matrix.
    memcpy(rC, iv, 16 * sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////
    
void Matrix4::transpose(void)
{
    double buf;
    buf = rC[1][0]; rC[1][0] = rC[0][1]; rC[0][1] = buf;
    buf = rC[2][0]; rC[2][0] = rC[0][2]; rC[0][2] = buf;
    buf = rC[3][0]; rC[3][0] = rC[0][3]; rC[0][3] = buf;
    buf = rC[2][1]; rC[2][1] = rC[1][2]; rC[1][2] = buf;
    buf = rC[3][1]; rC[3][1] = rC[1][3]; rC[1][3] = buf;
    buf = rC[3][2]; rC[3][2] = rC[2][3]; rC[2][3] = buf;
}

////////////////////////////////////////////////////////////////////////////////

AABBox3 smath::operator* (Matrix4 const & m, AABBox3 const & bb)
{
    // Create corner points.
    Point3 p1; p1.set(bb.getXMin(), bb.getYMin(), bb.getZMin());
    Point3 p2; p2.set(bb.getXMin(), bb.getYMin(), bb.getZMax());
    Point3 p3; p3.set(bb.getXMax(), bb.getYMin(), bb.getZMax());
    Point3 p4; p4.set(bb.getXMax(), bb.getYMin(), bb.getZMin());
    Point3 p5; p5.set(bb.getXMin(), bb.getYMax(), bb.getZMin());
    Point3 p6; p6.set(bb.getXMin(), bb.getYMax(), bb.getZMax());
    Point3 p7; p7.set(bb.getXMax(), bb.getYMax(), bb.getZMax());
    Point3 p8; p8.set(bb.getXMax(), bb.getYMax(), bb.getZMin());
    
    // Transform corner points.
    p1 = m * p1;
    p2 = m * p2;
    p3 = m * p3;
    p4 = m * p4;
    p5 = m * p5;
    p6 = m * p6;
    p7 = m * p7;
    p8 = m * p8;
    
    // Create transformed boundary box.
    AABBox3 bb2;
    bb2.grow(p1);
    bb2.grow(p2);
    bb2.grow(p3);
    bb2.grow(p4);
    bb2.grow(p5);
    bb2.grow(p6);
    bb2.grow(p7);
    bb2.grow(p8);
    return bb2;
}

////////////////////////////////////////////////////////////////////////////////

Matrix4 smath::operator* (Matrix4 const & m1, Matrix4 const & m2)
{
    Matrix4 m;
    
    m.rC[0][0] = (m1.rC[0][0] * m2.rC[0][0]) +
                 (m1.rC[0][1] * m2.rC[1][0]) +
                 (m1.rC[0][2] * m2.rC[2][0]) +
                 (m1.rC[0][3] * m2.rC[3][0]);
    m.rC[0][1] = (m1.rC[0][0] * m2.rC[0][1]) +
                 (m1.rC[0][1] * m2.rC[1][1]) +
                 (m1.rC[0][2] * m2.rC[2][1]) +
                 (m1.rC[0][3] * m2.rC[3][1]);
    m.rC[0][2] = (m1.rC[0][0] * m2.rC[0][2]) +
                 (m1.rC[0][1] * m2.rC[1][2]) +
                 (m1.rC[0][2] * m2.rC[2][2]) +
                 (m1.rC[0][3] * m2.rC[3][2]);
    m.rC[0][3] = (m1.rC[0][0] * m2.rC[0][3]) +
                 (m1.rC[0][1] * m2.rC[1][3]) +
                 (m1.rC[0][2] * m2.rC[2][3]) +
                 (m1.rC[0][3] * m2.rC[3][3]);
                 
    m.rC[1][0] = (m1.rC[1][0] * m2.rC[0][0]) +
                 (m1.rC[1][1] * m2.rC[1][0]) +
                 (m1.rC[1][2] * m2.rC[2][0]) +
                 (m1.rC[1][3] * m2.rC[3][0]);
    m.rC[1][1] = (m1.rC[1][0] * m2.rC[0][1]) +
                 (m1.rC[1][1] * m2.rC[1][1]) +
                 (m1.rC[1][2] * m2.rC[2][1]) +
                 (m1.rC[1][3] * m2.rC[3][1]);
    m.rC[1][2] = (m1.rC[1][0] * m2.rC[0][2]) +
                 (m1.rC[1][1] * m2.rC[1][2]) +
                 (m1.rC[1][2] * m2.rC[2][2]) +
                 (m1.rC[1][3] * m2.rC[3][2]);
    m.rC[1][3] = (m1.rC[1][0] * m2.rC[0][3]) +
                 (m1.rC[1][1] * m2.rC[1][3]) +
                 (m1.rC[1][2] * m2.rC[2][3]) +
                 (m1.rC[1][3] * m2.rC[3][3]);
                 
    m.rC[2][0] = (m1.rC[2][0] * m2.rC[0][0]) +
                 (m1.rC[2][1] * m2.rC[1][0]) +
                 (m1.rC[2][2] * m2.rC[2][0]) +
                 (m1.rC[2][3] * m2.rC[3][0]);
    m.rC[2][1] = (m1.rC[2][0] * m2.rC[0][1]) +
                 (m1.rC[2][1] * m2.rC[1][1]) +
                 (m1.rC[2][2] * m2.rC[2][1]) +
                 (m1.rC[2][3] * m2.rC[3][1]);
    m.rC[2][2] = (m1.rC[2][0] * m2.rC[0][2]) +
                 (m1.rC[2][1] * m2.rC[1][2]) +
                 (m1.rC[2][2] * m2.rC[2][2]) +
                 (m1.rC[2][3] * m2.rC[3][2]);
    m.rC[2][3] = (m1.rC[2][0] * m2.rC[0][3]) +
                 (m1.rC[2][1] * m2.rC[1][3]) +
                 (m1.rC[2][2] * m2.rC[2][3]) +
                 (m1.rC[2][3] * m2.rC[3][3]);
                 
    m.rC[3][0] = (m1.rC[3][0] * m2.rC[0][0]) +
                 (m1.rC[3][1] * m2.rC[1][0]) +
                 (m1.rC[3][2] * m2.rC[2][0]) +
                 (m1.rC[3][3] * m2.rC[3][0]);
    m.rC[3][1] = (m1.rC[3][0] * m2.rC[0][1]) +
                 (m1.rC[3][1] * m2.rC[1][1]) +
                 (m1.rC[3][2] * m2.rC[2][1]) +
                 (m1.rC[3][3] * m2.rC[3][1]);
    m.rC[3][2] = (m1.rC[3][0] * m2.rC[0][2]) +
                 (m1.rC[3][1] * m2.rC[1][2]) +
                 (m1.rC[3][2] * m2.rC[2][2]) +
                 (m1.rC[3][3] * m2.rC[3][2]);
    m.rC[3][3] = (m1.rC[3][0] * m2.rC[0][3]) +
                 (m1.rC[3][1] * m2.rC[1][3]) +
                 (m1.rC[3][2] * m2.rC[2][3]) +
                 (m1.rC[3][3] * m2.rC[3][3]);
                 
    return m;
}

////////////////////////////////////////////////////////////////////////////////
    
Normal3 smath::operator* (Matrix4 const & m, Normal3 const & n)
{
    // Compute inverse matrix.
    Matrix4 itm(m);
    itm.invert();
    itm.transpose();
    Normal3 r((itm.rC[0][0] * n.getX()) + (itm.rC[0][1] * n.getY()) +
              (itm.rC[0][2] * n.getZ()),
              (itm.rC[1][0] * n.getX()) + (itm.rC[1][1] * n.getY()) +
              (itm.rC[1][2] * n.getZ()),
              (itm.rC[2][0] * n.getX()) + (itm.rC[2][1] * n.getY()) +
              (itm.rC[2][2] * n.getZ()));
    return r;
}
        
////////////////////////////////////////////////////////////////////////////////
    
Point3 smath::operator* (Matrix4 const & m, Point3 const & p)
{
    Point3 r;
    double f = (m.rC[3][0] * p.getX()) + (m.rC[3][1] * p.getY()) +
               (m.rC[3][2] * p.getZ()) + (m.rC[3][3]);
    
    assert(smath::isSmallerEps(f) == false);
    
    r.setX(((m.rC[0][0] * p.getX()) + (m.rC[0][1] * p.getY()) +
            (m.rC[0][2] * p.getZ()) + (m.rC[0][3])) / f);
    r.setY(((m.rC[1][0] * p.getX()) + (m.rC[1][1] * p.getY()) +
            (m.rC[1][2] * p.getZ()) + (m.rC[1][3])) / f);
    r.setZ(((m.rC[2][0] * p.getX()) + (m.rC[2][1] * p.getY()) +
            (m.rC[2][2] * p.getZ()) + (m.rC[2][3])) / f);
        
    return r;
}
    
////////////////////////////////////////////////////////////////////////////////
    
Vector3 smath::operator* (Matrix4 const & m, Vector3 const & v)
{
    Vector3 r;
    
    r.setX((m.rC[0][0] * v.getX()) + (m.rC[0][1] * v.getY()) +
           (m.rC[0][2] * v.getZ()));
    r.setY((m.rC[1][0] * v.getX()) + (m.rC[1][1] * v.getY()) +
           (m.rC[1][2] * v.getZ()));
    r.setZ((m.rC[2][0] * v.getX()) + (m.rC[2][1] * v.getY()) +
           (m.rC[2][2] * v.getZ()));
    
    return r;
}

////////////////////////////////////////////////////////////////////////////////

std::istream & smath::operator>> (std::istream & is, Matrix4 & m)
{
    // Buffer to store the components; m should remain untouched if the
    // read operation fails!
    double n[16];
    // Initialize c; if >> fails, c should not be set!
    char c = 0;
    
    // Read (skipping whitespace).
    is >> c;
    if (c != '(') goto error;
    for (int i = 0; i < 15; ++i)
    {
        // Read number and separator.
        is >> n[i] >> c;
        // Expect comma.
        if (c != ',') goto error;
    }
    is >> n[15] >> c;
    if (c != ')') goto error;
    
    // Construct the matrix.
    m = Matrix4(n);
    
final:
    return is;
    
error:
    is.clear(std::ios_base::badbit);
    goto final;
}

////////////////////////////////////////////////////////////////////////////////

std::ostream & smath::operator<< (std::ostream & os, Matrix4 const & m)
{
    os << "("  << m.rC[0][0] // Row 1
       << ", " << m.rC[0][1]
       << ", " << m.rC[0][2]
       << ", " << m.rC[0][3]
       << ", " << m.rC[1][0] // Row 2
       << ", " << m.rC[1][1]
       << ", " << m.rC[1][2]
       << ", " << m.rC[1][3]
       << ", " << m.rC[2][0] // Row 3
       << ", " << m.rC[2][1]
       << ", " << m.rC[2][2]
       << ", " << m.rC[2][3]
       << ", " << m.rC[3][0] // Row 4
       << ", " << m.rC[3][1]
       << ", " << m.rC[3][2]
       << ", " << m.rC[3][3]
       << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////

// END
