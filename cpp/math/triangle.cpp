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

// STL headers.
#include <cassert>
#include <cmath>

// STEPS headers.
#include "../common.h"
#include "triangle.hpp"

////////////////////////////////////////////////////////////////////////////////

static inline void cross(double * v0, double * v1, double * vo)
{
    vo[0] = (v0[1] * v1[2]) - (v0[2] * v1[1]);
    vo[1] = (v0[2] * v1[0]) - (v0[0] * v1[2]);
    vo[2] = (v0[0] * v1[1]) - (v0[1] * v1[0]);
}

////////////////////////////////////////////////////////////////////////////////

double steps::math::triArea(double * v0, double * v1, double * v2)
{
    double vv[3];
    vv[0] = v1[0] - v0[0];
    vv[1] = v1[1] - v0[1];
    vv[2] = v1[2] - v0[2];
    double ww[3];
    ww[0] = v2[0] - v0[0];
    ww[1] = v2[1] - v0[1];
    ww[2] = v2[2] - v0[2];
    double cc[3];
    cross(vv, ww, cc);
    return 0.5 * sqrt((cc[0] * cc[0]) + (cc[1] * cc[1]) + (cc[2] * cc[2]));
}

////////////////////////////////////////////////////////////////////////////////

void steps::math::triBarycenter
(
    double * v0, double * v1, double * v2,
    double * po
)
{
    po[0] = (v0[0] + v1[0] + v2[0]) / 3.0;
    po[1] = (v0[1] + v1[1] + v2[1]) / 3.0;
    po[2] = (v0[2] + v1[2] + v2[2]) / 3.0;
}

////////////////////////////////////////////////////////////////////////////////

void steps::math::triNormal
(
    double * v0, double * v1, double * v2,
    double * vo
)
{
    double vv[3];
    vv[0] = v1[0] - v0[0];
    vv[1] = v1[1] - v0[1];
    vv[2] = v1[2] - v0[2];
    double ww[3];
    ww[0] = v2[0] - v0[0];
    ww[1] = v2[1] - v0[1];
    ww[2] = v2[2] - v0[2];
    double cc[3];
    cross(vv, ww, cc);
    double norm = sqrt((cc[0] * cc[0]) + (cc[1] * cc[1]) + (cc[2] * cc[2]));
    vo[0] = cc[0] / norm;
    vo[1] = cc[1] / norm;
    vo[2] = cc[2] / norm;
}

////////////////////////////////////////////////////////////////////////////////

// END
