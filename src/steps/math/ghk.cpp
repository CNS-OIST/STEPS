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

#include "ghk.hpp"

#include <cassert>
#include <cmath>

#include <easylogging++.h>

#include "constants.hpp"

#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////
/*
// The below may be applicable to chord conductance, but we assume a slope conductance
double steps::math::permeability
(
    double G, double V, int z, double T, double iconc, double oconc
)
{
    AssertLog(G >= 0.0);
    AssertLog(z != 0);
    AssertLog(T >= 0.0);
    AssertLog(iconc >= 0.0);
    AssertLog(oconc >= 0.0);

    double numerator = G * GAS_CONSTANT * T * (1-exp(((0-z)*V*FARADAY)/(GAS_CONSTANT*T)));
    double denominator = pow(z, 2.0) * pow(FARADAY, 2.0) * (iconc - (oconc*exp(((0-z)*V*FARADAY)/(GAS_CONSTANT*T))));

    return (numerator/denominator);
}
*/

double steps::math::permeability
(
    double G, double V, int z, double T, double iconc, double oconc
)
{
    AssertLog(G >= 0.0);
    AssertLog(z != 0);
    AssertLog(T >= 0.0);
    AssertLog(iconc >= 0.0);
    AssertLog(oconc >= 0.0);

    double B = iconc;
    double C = oconc;
    double D = ((0.0-static_cast<double>(z))*FARADAY)/(GAS_CONSTANT*T);

    double denominator = C*exp(2.0*D*V) + ((B-C)*D*V -C -B)*exp(D*V) + B;
    double numerator = exp(2.0*D*V) - 2.0*exp(D*V) + 1;

    double A = G * (numerator/denominator);

    double P  = (A*GAS_CONSTANT*T)/(pow(z, 2.0)*pow(FARADAY, 2.0));

    return P;

}

////////////////////////////////////////////////////////////////////////////////

double steps::math::GHKcurrent
(
    double P, double V, int z, double T, double iconc, double oconc
)
{
    AssertLog(z != 0);
    AssertLog(T >= 0.0);
    AssertLog(iconc >= 0.0);
    AssertLog(oconc >= 0.0);

    double numerator = (P * pow(z, 2.0) * V * pow(FARADAY, 2.0))/(GAS_CONSTANT*T) * (iconc - (oconc*exp(((0-z)*V*FARADAY)/(GAS_CONSTANT*T))));
    double denominator = (1-exp(((0-z)*V*FARADAY)/(GAS_CONSTANT*T)));

    return (numerator/denominator);
}
