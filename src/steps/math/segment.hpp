/*
 ___license_placeholder___
 */

/** \file math/segment.hpp
 *  Geometric functions on segments in 3-d space.
 */

#pragma once

#include "point.hpp"

namespace steps::math {

/** Intersects a line segment with a point
 *
 * \param a, b: Vertices of the ray/segment.
 * \param p: Point.
 * \param is_segment: When true, the two points delimit a segment.
 *     Otherwise they define a ray starting at a and passing by b.
 * \return True if the point lies on the ray/segment. False otherwise.
 */
bool segment_intersect_point(const point3d& a,
                             const point3d& b,
                             const point3d& p,
                             bool is_segment = true);

/** Intersects a line segment with a line segment
 *
 * There is an assert to prevent colinearity. In general
 * the 2 segments should not be colinear.
 *
 * \param p0, p1: Vertices of the ray/segment p.
 * \param q0, q1: Vertices of the ray/segment q.
 * \param intersection0: Intersection point and possible start of intersection line.
 * \param intersection1: Intersection point and possible end of intersection line.
 * \param is_segment_p: When true, the two points delimit a segment.
 *     Otherwise they define a ray starting at p0 and passing by p1.
 * \param is_segment_q: When true, the two points delimit a segment.
 *     Otherwise they define a ray starting at q0 and passing by q1.
 *
 * \return True if they intersect in a single point.
 */
bool segment_intersect_segment(const point3d& p0,
                               const point3d& p1,
                               const point3d& q0,
                               const point3d& q1,
                               point3d& intersection0,
                               point3d& intersection1,
                               const bool is_segment_p = true,
                               const bool is_segment_q = true);

}  // namespace steps::math
