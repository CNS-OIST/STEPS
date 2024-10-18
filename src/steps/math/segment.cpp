/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#include "segment.hpp"
#include "triangle.hpp"
#include "util/error.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace steps::math {

bool segment_intersect_point(const point3d& a,
                             const point3d& b,
                             const point3d& p,
                             const bool is_segment) {
    // Compute if point lies of a segment paying attention to numerical precision problems

    // 4 ULP
    constexpr auto adimensional_tol = 4 * std::numeric_limits<point3d::value_type>::epsilon();
    const auto tol = adimensional_tol * normL1(a - b);

    if (a.almostEqual(p, tol) || b.almostEqual(p, tol)) {
        return true;
    }
    if (a.almostEqual(b, tol)) {
        return false;
    }

    const point3d ab = b - a;
    const point3d ap = p - a;
    const double denom = ab.dot(ab);
    const double toldenom = adimensional_tol * denom;
    const double tdenom = ab.dot(ap);

    if (is_segment && (tdenom < -toldenom || tdenom > toldenom + denom)) {
        return false;
    }

    // projection on ab
    const point3d q = a + (tdenom / denom) * ab;
    return p.almostEqual(q, tol);
}

bool segment_intersect_segment(const point3d& p0,
                               const point3d& p1,
                               const point3d& q0,
                               const point3d& q1,
                               point3d& intersection0,
                               point3d& intersection1,
                               const bool is_segment_p,
                               const bool is_segment_q) {
    const auto u = p1 - p0;
    const auto v = q1 - q0;

    // 4 ULP
    constexpr auto adimensional_tol = 4 * std::numeric_limits<point3d::value_type>::epsilon();
    const auto max_l = std::max(normL1(u), normL1(v));
    const auto tol = adimensional_tol * max_l;

    // p0 - p1 is not a real segment
    if (p0.almostEqual(p1, tol)) {
        if (segment_intersect_point(q0, q1, p0)) {
            intersection0 = p0;
            intersection1 = intersection0;
            return true;
        }
        return false;
    }

    // q0 - q1 is not a real segment
    if (q0.almostEqual(q1, tol)) {
        if (segment_intersect_point(p0, p1, q0)) {
            intersection0 = q0;
            intersection1 = intersection0;
            return true;
        }
        return false;
    }

    const auto w = p0 - q0;

    const double a = squared_norm(u);
    const double b = u.dot(v);
    const double c = squared_norm(v);
    const double d = w.dot(u);
    const double e = w.dot(v);

    const double denom = a * c - b * b;
    if (std::abs(denom) < adimensional_tol * max_l * max_l * max_l * max_l) {
        // check that they do not overlap
        if (norm(v.cross(w)) > (tol * std::sqrt(c))) {
            return false;
        }

        // they overlap find intersection line
        if (!is_segment_p && !is_segment_q) {
            intersection0 = point3d{};
            intersection1 = point3d{};
            CLOG(WARNING, "general_log")
                << "Intersection of colinear lines. Intersection points are meaningless.\n";
            return true;
        } else if (!is_segment_p) {
            intersection0 = q0;
            intersection1 = q1;
            return true;
        } else if (!is_segment_q) {
            intersection0 = p0;
            intersection1 = p1;
            return true;
        } else {
            // find q0 and q1 position based on on p0 -p1 coordinates
            const auto z = q1 - p0;
            // w is defined inverted compared to the other vectors
            const double tq0 = -w.dot(u) / a;
            const double tq1 = z.dot(u) / a;
            std::vector<std::pair<double, point3d>> stor{{tq0, q0},
                                                         {tq1, q1},
                                                         {0.0, p0},
                                                         {1.0, p1}};

            // make stor[0] the smallest one. Always
            if (tq0 > tq1) {
                std::swap(stor[0], stor[1]);
            }

            // along the line, after p1
            if (stor[0].first > 1.0) {
                return false;
            }

            // along the line, before p0
            if (stor[1].first < 0.0) {
                return false;
            }

            // one after the other
            if (stor[0].second.almostEqual(p1, tol)) {
                intersection0 = p1;
                intersection1 = intersection0;
                return true;
            }

            // one after the other
            if (stor[1].second.almostEqual(p0, tol)) {
                intersection0 = p0;
                intersection1 = intersection0;
                return true;
            }

            std::sort(stor.begin(), stor.end(), [](const auto& val0, const auto& val1) {
                return val0.first < val1.first;
            });
            // we take always the internal intersection
            intersection0 = stor[1].second;
            intersection1 = stor[2].second;
            return true;
        }
    }
    const double denomtol = denom * adimensional_tol;

    // handle extremes
    if (segment_intersect_point(p0, p1, q0)) {
        intersection0 = q0;
        intersection1 = intersection0;
        return true;
    }
    if (segment_intersect_point(p0, p1, q1)) {
        intersection0 = q1;
        intersection1 = intersection0;
        return true;
    }
    if (segment_intersect_point(q0, q1, p0)) {
        intersection0 = p0;
        intersection1 = intersection0;
        return true;
    }
    if (segment_intersect_point(q0, q1, p1)) {
        intersection0 = p1;
        intersection1 = intersection0;
        return true;
    }

    const double tdenom = (b * e - c * d);
    if (is_segment_p && (tdenom < -denomtol || tdenom > denom + denomtol)) {
        return false;
    }
    const double sdenom = (a * e - b * d);
    if (is_segment_q && (sdenom < -denomtol || sdenom > denom + denomtol)) {
        return false;
    }

    // if on the same plane, the segments intersect
    if (tri_intersect_point(p0, p1, q0, q1, false)) {
        intersection0 = (p0 + (tdenom / denom) * u + q0 + (sdenom / denom) * v) / 2.0;
        intersection1 = intersection0;
        return true;
    }

    return false;
}

}  // namespace steps::math
