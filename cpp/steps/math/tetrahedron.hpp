////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_MATH_TETRAHEDRON_HPP
#define STEPS_MATH_TETRAHEDRON_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>

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
