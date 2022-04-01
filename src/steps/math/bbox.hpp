/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

#pragma once

#include <utility>

#include "point.hpp"

namespace steps {
namespace math {

/** Axis-aligned closed bounding box.
 */

struct bounding_box {
    /** Construct empty bounding box. */
    bounding_box(): empty_flag(true) {}

    /** Construct bounding box encompassing a single point. */
    explicit bounding_box(const point3d &x):
        x0(x), x1(x), empty_flag(false) {}

    /** Construct bounding box representing the interval [x0,x1]
     * with respect to the partial order on points. */
    bounding_box(const point3d &x0_, const point3d &x1_):
        x0(x0_), x1(x1_), empty_flag(!partial_leq(x0_,x1_)) {}

    /** Construct bounding box that encompasses the points descibed
     * by given iterator range. */
    template <typename I>
    bounding_box(I b, I e): empty_flag(true) {
        while (b != e) insert(*b++);
    }

    /** Expand bounding box to minimally enclose point p. */
    void insert(const point3d &p) {
        if (empty()) {
            x0 = x1 = p;
            empty_flag = false;
        }
        else {
            x0 = meet(x0, p);
            x1 = join(x1, p);
        }
    }

    /* Reset bounding box to empty state. */
    void clear() { empty_flag = true; }

    /* Return true if bounding box is empty. */
    bool empty() const { return empty_flag; }

    /* Return true if point p lies within bounding box. */
    bool contains(const point3d &p) const {
        return !empty() && partial_leq(x0,p) && partial_leq(p,x1);
    }

    /** Return reference to minimum point of bounding box.
     *
     * Note that value not well-defined if bounding box is empty.
     */
    const point3d &min() const { return x0; }

    /** Return reference to maximum point of bounding box.
     *
     * Note that value not well-defined if bounding box is empty.
     */
    const point3d &max() const { return x1; }

private:
    point3d x0, x1;
    bool empty_flag;

    static point3d meet(const point3d &a, const point3d &b) {
        return point3d{std::min(a[0],b[0]), std::min(a[1],b[1]), std::min(a[2],b[2])};
    }

    static point3d join(const point3d &a, const point3d &b) {
        return point3d{std::max(a[0],b[0]), std::max(a[1],b[1]), std::max(a[2],b[2])};
    }

    static bool partial_leq(const point3d &a, const point3d &b) {
        return a[0]<=b[0] && a[1]<=b[1] && a[2]<=b[2];
    }
};


} // namespace math
} // namespace steps
