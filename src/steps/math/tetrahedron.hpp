/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

/** /file Geometric functions on tetrahedra in 3-d space.
 */

#ifndef STEPS_MATH_TETRAHEDRON_HPP
#define STEPS_MATH_TETRAHEDRON_HPP 1

// STEPS headers.
#include "steps/common.h"
#include "steps/math/point.hpp"

namespace steps {
namespace math {

/** Calculate volume of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Volume.
 */
double tet_vol(const point3d &p0, const point3d &p1,
               const point3d &p2, const point3d &p3);

/** Calculate tetrahedron barycenter.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Barycenter.
 */
point3d tet_barycenter(const point3d &p0, const point3d &p1,
                       const point3d &p2, const point3d &p3);

/** Test for point inclusion in tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \param pi Point to test
 * \return True if pi lies within tetrahedron.
 */
bool tet_inside(const point3d &p0, const point3d &p1,
                const point3d &p2, const point3d &p3,
                const point3d &pi);

/** Select point in tetrahedron from uniformly generated variates.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \param s,t,u Three uniformly-generated variates in [0,1].
 * \return Sampled point.
 */
point3d tet_ranpnt(const point3d &p0, const point3d &p1,
                   const point3d &p2, const point3d &p3,
                   double s, double t, double u);

/** Caclulate circumradius of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Circumradius.
 */
double tet_circumrad(const point3d &p0, const point3d &p1,
                     const point3d &p2, const point3d &p3);

/** Caclulate square of circumradius of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Circumradius squared.
 */
double tet_circumrad2(const point3d &p0, const point3d &p1,
                      const point3d &p2, const point3d &p3);

/** Caclulate length of shortest edge of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Shortest edge length.
 */
double tet_shortestedge(const point3d &p0, const point3d &p1,
                        const point3d &p2, const point3d &p3);

/** Caclulate squared length of shortest edge of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Shortest edge length squared.
 */
double tet_shortestedge2(const point3d &p0, const point3d &p1,
                         const point3d &p2, const point3d &p3);

/** Caclulate length of longest edge of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Shortest edge length.
 */
double tet_longestedge(const point3d &p0, const point3d &p1,
                        const point3d &p2, const point3d &p3);

/** Caclulate squared length of longest edge of tetrahedron.
 *
 * \param p0,p1,p2,p3 Vertices of tetrahedron.
 * \return Shortest edge length squared.
 */
double tet_longestedge2(const point3d &p0, const point3d &p1,
                         const point3d &p2, const point3d &p3);

}} // namespace steps:math

#endif // ndef STEPS_MATH_TETRAHEDRON_HPP
// END
