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

#ifndef STEPS_MATH_TETRAHEDRON_HPP
#define STEPS_MATH_TETRAHEDRON_HPP 1


// STEPS headers.
#include "../common.h"

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN double tet_vol
(
    double * v0, double * v1,
    double * v2, double * v3
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN void tet_barycenter
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * po
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN void tet_barycentric
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * pi, double * po
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN bool tet_inside
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * pi
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN void tet_ranpnt
(
    double * v0, double * v1,
    double * v2, double * v3,
    double r_unf0, double r_unf1, double r_unf2,
    double * po
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN void tet_circumsphere
(
    double * v0, double * v1,
    double * v2, double * v3,
    double * ctr, double * rad
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN double tet_circumrad
(
    double * v0, double * v1,
    double * v2, double * v3
);

STEPS_EXTERN double tet_circumrad2
(
    double * v0, double * v1,
    double * v2, double * v3
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN double tet_shortestedge
(
    double * v0, double * v1,
    double * v2, double * v3
);

STEPS_EXTERN double tet_shortestedge2
(
    double * v0, double * v1,
    double * v2, double * v3
);

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN double tet_longestedge
(
    double * v0, double * v1,
    double * v2, double * v3
);

STEPS_EXTERN double tet_longestedge2
(
    double * v0, double * v1,
    double * v2, double * v3
);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_TETRAHEDRON_HPP

// END
