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

#ifndef STEPS_MATH_CONSTANTS_HPP
#define STEPS_MATH_CONSTANTS_HPP 1


// STEPS headers.
#include "../common.h"

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

const double E                          = 2.71828182845904523536028747135;
const double PI                         = 3.14159265358979323846264338328;

const float  IEEE_MACH_EPSILON32        = 2e-24;
const double IEEE_MACH_EPSILON64        = 2e-53;
const float  IEEE_EPSILON32             = 1.0e-7;
const double IEEE_EPSILON64             = 1.0e-15;

const double  IEEE_HUGE                 = 1e150;

////////////////////////////////////////////////////////////////////////////////

// Source: physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadrp
const double AVOGADRO                   = 6.0221415e23;

////////////////////////////////////////////////////////////////////////////////

const double M_TO_NM                    = 1.0e9;
const double M_TO_UM                    = 1.0e6;
const double M_TO_MM                    = 1.0e3;
const double M_TO_CM                    = 1.0e2;
const double M_TO_DM                    = 1.0e1;

const double UM_TO_M                    = 1.0e-6;

////////////////////////////////////////////////////////////////////////////////

const double M2_TO_UM2                  = 1.0e12;

const double UM2_TO_M2                  = 1.0e-12;

////////////////////////////////////////////////////////////////////////////////

const double M3_TO_DM3                  = 1.0e3;

const double DM3_TO_M3                  = 1.0e-3;

const double M3_TO_UM3                  = 1.0e18;

const double UM3_TO_M3                  = 1.0e-18;

////////////////////////////////////////////////////////////////////////////////

const double MS_TO_S                    = 1.0e-3;

const double S_TO_MS                    = 1.0e3;

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_CONSTANTS_HPP

// END
