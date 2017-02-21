/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_MATH_GHK_HPP
#define STEPS_MATH_GHK_HPP 1

// STEPS headers.
#include "steps/common.h"

 namespace steps {
 namespace math {

////////////////////////////////////////////////////////////////////////////////

// Return the permeability in the GHK flux equation from given values of:
// G (slope conductance in siemens), V (voltage in volts), z (valence),
// T (temperature in kelvin),
// iconc (inner concentration of ion in mol per cubic meter),
// oconc (outer concentration of ion in mol per cubic meter)

STEPS_EXTERN double permeability
(
    double G, double V, int z, double T, double iconc, double oconc
);

////////////////////////////////////////////////////////////////////////////////

// Return the single-channel current from the GHK flux equation from given
// value of:
// P (single-channel permeability in meters cubed/second), V (voltage in volts),
// z (valence), T (temperature in kelvin),
// iconc (inner concentration of ion in mol per cubic meter),
// oconc (outer concentration of ion in mol per cubic meter)

STEPS_EXTERN double GHKcurrent
(
    double P, double V, int z, double T, double iconc, double oconc
);

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MATH_GHK_HPP

// END
