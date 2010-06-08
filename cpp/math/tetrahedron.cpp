////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

// STEPS headers.
#include "../common.h"
#include "linsolve.hpp"
#include "tetrahedron.hpp"

// STL headers.
#include <cassert>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////

static inline double det4X3
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    return v1[2]*v2[1]*v3[0]-v0[2]*v2[1]*v3[0]-1.0*v1[1]*v2[2]*v3[0]+v0[1]*
           v2[2]*v3[0]+v0[2]*v1[1]*v3[0]-v0[1]*v1[2]*v3[0]-1.0*v1[2]*v2[0]*
           v3[1]+v0[2]*v2[0]*v3[1]+1.0*v1[0]*v2[2]*v3[1]-v0[0]*v2[2]*
           v3[1]-v0[2]*v1[0]*v3[1]+v0[0]*v1[2]*v3[1]+1.0*v1[1]*v2[0]*
           v3[2]-v0[1]*v2[0]*v3[2]-1.0*v1[0]*v2[1]*v3[2]+v0[0]*v2[1]*
           v3[2]+v0[1]*v1[0]*v3[2]-v0[0]*v1[1]*v3[2]-v0[2]*v1[1]*v2[0]*
           1.0+v0[1]*v1[2]*v2[0]*1.0+v0[2]*v1[0]*v2[1]*1.0-v0[0]*v1[2]*
           v2[1]*1.0-v0[1]*v1[0]*v2[2]*1.0+v0[0]*v1[1]*v2[2];
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_vol
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    return fabs(det4X3(v0, v1, v2, v3) / 6.0);
}

////////////////////////////////////////////////////////////////////////////////

void steps::math::tet_barycenter
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * po
)
{
    po[0] = (v0[0] + v1[0] + v2[0] + v3[0]) / 4.0;
    po[1] = (v0[1] + v1[1] + v2[1] + v3[1]) / 4.0;
    po[2] = (v0[2] + v1[2] + v2[2] + v3[2]) / 4.0;
}

////////////////////////////////////////////////////////////////////////////////

void steps::math::tet_barycentric
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * pi, double * po
)
{
    // Copy first tet corner point in local variables.
    double x0 = v0[0];
    double y0 = v0[1];
    double z0 = v0[2];
    // Prepase base matrix.
    double a[12];
    a[ 0] = v1[0] - x0;
    a[ 1] = v1[1] - y0;
    a[ 2] = v1[2] - z0;
    a[ 3] = v2[0] - x0;
    a[ 4] = v2[1] - y0;
    a[ 5] = v2[2] - z0;
    a[ 6] = v3[0] - x0;
    a[ 7] = v3[1] - y0;
    a[ 8] = v3[2] - z0;
    a[ 9] = pi[0] - x0;
    a[10] = pi[1] - y0;
    a[11] = pi[2] - z0;

    // Solve this system.
    int info = steps::math::linsolve(3, 1, a);
    if (info != 0)
    {
        // TODO: error message
        assert(0);
    }

    // Build results.
    double b1 = a[ 9];
    double b2 = a[10];
    double b3 = a[11];
    if (b2 < b1)
    {
        double temp = b2;
        b2 = b1;
        b1 = temp;
    }
    if (b3 < b2)
    {
        double temp = b3;
        b3 = b2;
        b2 = temp;
        if (b2 < b1)
        {
            double temp = b2;
            b2 = b1;
            b1 = temp;
        }
    }
    po[0] = 1.0 - (b1 + b2 + b3);
    po[1] = b1;
    po[2] = b2;
    po[3] = b3;
}

////////////////////////////////////////////////////////////////////////////////

bool steps::math::tet_inside
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * pi
)
{
    double po[4];
    steps::math::tet_barycentric(v0, v1, v2, v3, pi, po);
    if (po[0] < 0.0) return false;
    if (po[1] < 0.0) return false;
    if (po[2] < 0.0) return false;
    if (po[3] < 0.0) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void steps::math::tet_ranpnt
(
    double * v0, double * v1,
    double * v2, double * v3,
    double s, double t, double u,
    double * po
)
{
    if (s + t > 1.0)
    {
        // cut'n fold the cube into a prism
        s = 1.0 - s;
        t = 1.0 - t;
    }
    if (t + u > 1.0)
    {
        // cut'n fold the prism into a tetrahedron
        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
    }
    else if (s + t + u > 1.0)
    {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
    }
    // a,s,t,u are the barycentric coordinates of the random point.
    double a = 1 - s - t - u;
    po[0] = (a * v0[0]) + (s * v1[0]) +
            (t * v2[0]) + (u * v3[0]);
    po[1] = (a * v0[1]) + (s * v1[1]) +
            (t * v2[1]) + (u * v3[1]);
    po[2] = (a * v0[2]) + (s * v1[2]) +
            (t * v2[2]) + (u * v3[2]);
}

////////////////////////////////////////////////////////////////////////////////

void steps::math::tet_circumsphere
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * ctr, double * rad
)
{
    double x0 = v0[0];
    double y0 = v0[1];
    double z0 = v0[2];

    double a[12];

    double pp = a[ 0] = v1[0] - x0;
    double qq = a[ 3] = v1[1] - y0;
    double ss = a[ 6] = v1[2] - z0;
    a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);

    pp = a[ 1] = v2[0] - x0;
    qq = a[ 4] = v2[1] - y0;
    ss = a[ 7] = v2[2] - z0;
    a[10] = (pp * pp) + (qq * qq) + (ss * ss);

    pp = a[ 2] = v3[0] - x0;
    qq = a[ 5] = v3[1] - y0;
    ss = a[ 8] = v3[2] - z0;
    a[11] = (pp * pp) + (qq * qq) + (ss * ss);

    int info = steps::math::linsolve(3, 1, a);
    if (info != 0)
    {
        assert(0);
    }

    double xx0 = a[ 9];
    double yy0 = a[10];
    double zz0 = a[11];

    ctr[0] = x0 + (0.5 * xx0);
    ctr[1] = y0 + (0.5 * yy0);
    ctr[2] = z0 + (0.5 * zz0);
    *rad = 0.5 * sqrt((xx0 * xx0) + (yy0 * yy0) + (zz0 * zz0));
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_circumrad
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    double x0 = v0[0];
    double y0 = v0[1];
    double z0 = v0[2];

    double a[12];

    double pp = a[ 0] = v1[0] - x0;
    double qq = a[ 3] = v1[1] - y0;
    double ss = a[ 6] = v1[2] - z0;
    a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);

    pp = a[ 1] = v2[0] - x0;
    qq = a[ 4] = v2[1] - y0;
    ss = a[ 7] = v2[2] - z0;
    a[10] = (pp * pp) + (qq * qq) + (ss * ss);

    pp = a[ 2] = v3[0] - x0;
    qq = a[ 5] = v3[1] - y0;
    ss = a[ 8] = v3[2] - z0;
    a[11] = (pp * pp) + (qq * qq) + (ss * ss);

    int info = steps::math::linsolve(3, 1, a);
    if (info != 0)
    {
        assert(0);
    }

    double xx0 = a[ 9];
    double yy0 = a[10];
    double zz0 = a[11];

    return 0.5 * sqrt((xx0 * xx0) + (yy0 * yy0) + (zz0 * zz0));
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_circumrad2
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    double x0 = v0[0];
    double y0 = v0[1];
    double z0 = v0[2];

    double a[12];

    double pp = a[ 0] = v1[0] - x0;
    double qq = a[ 3] = v1[1] - y0;
    double ss = a[ 6] = v1[2] - z0;
    a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);

    pp = a[ 1] = v2[0] - x0;
    qq = a[ 4] = v2[1] - y0;
    ss = a[ 7] = v2[2] - z0;
    a[10] = (pp * pp) + (qq * qq) + (ss * ss);

    pp = a[ 2] = v3[0] - x0;
    qq = a[ 5] = v3[1] - y0;
    ss = a[ 8] = v3[2] - z0;
    a[11] = (pp * pp) + (qq * qq) + (ss * ss);

    int info = steps::math::linsolve(3, 1, a);
    if (info != 0)
    {
        assert(0);
    }

    double xx0 = a[ 9];
    double yy0 = a[10];
    double zz0 = a[11];

    return 0.25 * ((xx0 * xx0) + (yy0 * yy0) + (zz0 * zz0));
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_shortestedge
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    // Edge 0: v1 - v0
    double dx = v1[0] - v0[0];
    double dy = v1[1] - v0[1];
    double dz = v1[2] - v0[2];
    double ls = (dx * dx) + (dy * dy) + (dz * dz);

    // Edge 1: v2 - v0
    dx = v2[0] - v0[0];
    dy = v2[1] - v0[1];
    dz = v2[2] - v0[2];
    double ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 2: v3 - v0
    dx = v3[0] - v0[0];
    dy = v3[1] - v0[1];
    dz = v3[2] - v0[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 3: v2 - v1
    dx = v2[0] - v1[0];
    dy = v2[1] - v1[1];
    dz = v2[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 4: v3 - v1
    dx = v3[0] - v1[0];
    dy = v3[1] - v1[1];
    dz = v3[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 5: v3 - v2
    dx = v3[0] - v2[0];
    dy = v3[1] - v2[1];
    dz = v3[2] - v2[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    return sqrt(ls);
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_shortestedge2
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    // Edge 0: v1 - v0
    double dx = v1[0] - v0[0];
    double dy = v1[1] - v0[1];
    double dz = v1[2] - v0[2];
    double ls = (dx * dx) + (dy * dy) + (dz * dz);

    // Edge 1: v2 - v0
    dx = v2[0] - v0[0];
    dy = v2[1] - v0[1];
    dz = v2[2] - v0[2];
    double ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 2: v3 - v0
    dx = v3[0] - v0[0];
    dy = v3[1] - v0[1];
    dz = v3[2] - v0[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 3: v2 - v1
    dx = v2[0] - v1[0];
    dy = v2[1] - v1[1];
    dz = v2[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 4: v3 - v1
    dx = v3[0] - v1[0];
    dy = v3[1] - v1[1];
    dz = v3[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    // Edge 5: v3 - v2
    dx = v3[0] - v2[0];
    dy = v3[1] - v2[1];
    dz = v3[2] - v2[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln < ls) ls = ln;

    return ls;
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_longestedge
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    // Edge 0: v1 - v0
    double dx = v1[0] - v0[0];
    double dy = v1[1] - v0[1];
    double dz = v1[2] - v0[2];
    double ll = (dx * dx) + (dy * dy) + (dz * dz);

    // Edge 1: v2 - v0
    dx = v2[0] - v0[0];
    dy = v2[1] - v0[1];
    dz = v2[2] - v0[2];
    double ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 2: v3 - v0
    dx = v3[0] - v0[0];
    dy = v3[1] - v0[1];
    dz = v3[2] - v0[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 3: v2 - v1
    dx = v2[0] - v1[0];
    dy = v2[1] - v1[1];
    dz = v2[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 4: v3 - v1
    dx = v3[0] - v1[0];
    dy = v3[1] - v1[1];
    dz = v3[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 5: v3 - v2
    dx = v3[0] - v2[0];
    dy = v3[1] - v2[1];
    dz = v3[2] - v2[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    return sqrt(ll);
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::tet_longestedge2
(
    double * v0, double * v1,
    double * v2, double * v3
)
{
    // Edge 0: v1 - v0
    double dx = v1[0] - v0[0];
    double dy = v1[1] - v0[1];
    double dz = v1[2] - v0[2];
    double ll = (dx * dx) + (dy * dy) + (dz * dz);

    // Edge 1: v2 - v0
    dx = v2[0] - v0[0];
    dy = v2[1] - v0[1];
    dz = v2[2] - v0[2];
    double ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 2: v3 - v0
    dx = v3[0] - v0[0];
    dy = v3[1] - v0[1];
    dz = v3[2] - v0[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 3: v2 - v1
    dx = v2[0] - v1[0];
    dy = v2[1] - v1[1];
    dz = v2[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 4: v3 - v1
    dx = v3[0] - v1[0];
    dy = v3[1] - v1[1];
    dz = v3[2] - v1[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    // Edge 5: v3 - v2
    dx = v3[0] - v2[0];
    dy = v3[1] - v2[1];
    dz = v3[2] - v2[2];
    ln = (dx * dx) + (dy * dy) + (dz * dz);
    if (ln > ll) ll = ln;

    return ll;
}

////////////////////////////////////////////////////////////////////////////////

// END

