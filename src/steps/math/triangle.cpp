/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

 */

#include <cmath>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/point.hpp"
#include "steps/math/triangle.hpp"

namespace steps {
namespace math {

////////////////////////////////////////////////////////////////////////////////

double tri_area(const point3d &p0, const point3d &p1, const point3d &p2)
{
    point3d c = cross(p1-p0, p2-p0);
    return 0.5 * std::sqrt(dot(c, c));
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_barycenter(const point3d &p0, const point3d &p1, const point3d &p2)
{
    return (p0+p1+p2)/3.0;
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_normal(const point3d &p0, const point3d &p1, const point3d &p2)
{
    point3d c = cross(p1-p0, p2-p0);
    return c/std::sqrt(dot(c,c));
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_ranpnt(const point3d &p0, const point3d &p1, const point3d &p2, double s, double t)
{
    // from http://www.cs.princeton.edu/~funk/tog02.pdf
    double u = std::sqrt(s);
    double v = u*t;

    return (1-u)*p0 + (u-v)*p1 + v*p2;
}

}} // namespace steps::math

// END
