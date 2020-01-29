/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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

/** /file Geometric functions on triangles in 3-d space.
 */

#ifndef STEPS_MATH_TRIANGLE_HPP
#define STEPS_MATH_TRIANGLE_HPP 1


// STEPS headers.
#include "steps/common.h"
#include "steps/math/point.hpp"

namespace steps {
namespace math {

/** Calculate area of triangle.
 *
 * \param p0,p1,p2 Vertices of triangle.
 * \return Area.
 */
double tri_area(const point3d &p0, const point3d &p1, const point3d &p2);

/** Calculate triangle barycentre.
 *
 * \param p0,p1,p2 Vertices of triangle.
 * \return Barycentre.
 */
point3d tri_barycenter(const point3d &p0, const point3d &p1, const point3d &p2);

/** Calculate triangle normal vector.
 *
 * \param p0,p1,p2 Vertices of triangle.
 * \return Unit length normal vector.
 */
point3d tri_normal(const point3d &p0, const point3d &p1, const point3d &p2);

/** Select point in triangle from uniformly generated variates.
 *
 * \param p0,p1,p2 Vertices of triangle.
 * \param s,t Two uniformly-generated variates in [0,1].
 * \return Sampled point.
 */
point3d tri_ranpnt(const point3d &p0, const point3d &p1, const point3d &p2, double s, double t);

/** Intersects a triangle with a line segment
 *
 * \param tp0, tp1, tp2: Vertices of triangle
 * \param lp0, lp1: Two points of the ray/segment
 * \param is_segment: When true, the two points delimit a segment.
 *     Otherwise they define a ray starting at lp0 and passing by lp1
 * \return The intersection point if exists, else
 */
bool tri_intersect_line(const point3d &tp0, const point3d &tp1, const point3d &tp2,
                        const point3d &lp0, const point3d &lp1, point3d &intersection,
                        bool is_segment=true);

}} // namespace steps::math

#endif // ndef STEPS_MATH_TRIANGLE_HPP
// END
