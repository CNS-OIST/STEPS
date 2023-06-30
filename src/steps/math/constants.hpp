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

#pragma once

#include "util/common.hpp"

namespace steps::math {

inline constexpr double E = 2.71828182845904523536028747135;
inline constexpr double PI = 3.14159265358979323846264338328;

inline constexpr float IEEE_MACH_EPSILON32 = 2e-24;
inline constexpr double IEEE_MACH_EPSILON64 = 2e-53;
inline constexpr float IEEE_EPSILON32 = 1.0e-7;
inline constexpr double IEEE_EPSILON64 = 1.0e-15;

inline constexpr double IEEE_HUGE = 1e150;

////////////////////////////////////////////////////////////////////////////////

// Avogadro constant
// Source: http://physics.nist.gov/cuu/Constants/index.html  4/3/2010
inline constexpr double AVOGADRO = 6.02214179e23;

// Faraday constant
// Source: http://physics.nist.gov/cuu/Constants/index.html     4/3/2010
inline constexpr double FARADAY = 96485.3399;

// Molar gas constant
// Source: http://physics.nist.gov/cuu/Constants/index.html  4/3/2010
inline constexpr double GAS_CONSTANT = 8.314472;

// The Elementray charge
// Source: http://physics.nist.gov/cuu/Constants/index.html  10/3/2010
inline constexpr double E_CHARGE = 1.602176487e-19;

////////////////////////////////////////////////////////////////////////////////

inline constexpr double M_TO_NM = 1.0e9;
inline constexpr double M_TO_UM = 1.0e6;
inline constexpr double M_TO_MM = 1.0e3;
inline constexpr double M_TO_CM = 1.0e2;
inline constexpr double M_TO_DM = 1.0e1;

inline constexpr double UM_TO_M = 1.0e-6;

////////////////////////////////////////////////////////////////////////////////

inline constexpr double M2_TO_UM2 = 1.0e12;

inline constexpr double UM2_TO_M2 = 1.0e-12;

////////////////////////////////////////////////////////////////////////////////

inline constexpr double M3_TO_DM3 = 1.0e3;

inline constexpr double DM3_TO_M3 = 1.0e-3;

inline constexpr double M3_TO_UM3 = 1.0e18;

inline constexpr double UM3_TO_M3 = 1.0e-18;

////////////////////////////////////////////////////////////////////////////////

inline constexpr double MS_TO_S = 1.0e-3;

inline constexpr double S_TO_MS = 1.0e3;

}  // namespace steps::math
