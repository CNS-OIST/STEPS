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

#include "ghk.hpp"

#include <cassert>
#include <cmath>


#include "constants.hpp"

#include "util/error.hpp"

namespace steps::math {

double permeability(double G, double V, int z, double T, double iconc, double oconc) {
    AssertLog(G >= 0.0);
    AssertLog(z != 0);
    AssertLog(T >= 0.0);
    AssertLog(iconc >= 0.0);
    AssertLog(oconc >= 0.0);

    double B = iconc * 1000.0;
    double C = oconc * 1000.0;
    double D = ((0.0 - static_cast<double>(z)) * FARADAY) / (GAS_CONSTANT * T);

    double denominator = C * exp(2.0 * D * V) + ((B - C) * D * V - C - B) * exp(D * V) + B;
    double numerator = exp(2.0 * D * V) - 2.0 * exp(D * V) + 1;

    // If both the denominator and numerator are 0, we compute the limit of the indeterminate form
    // by applying L’Hôpital’s rule:
    double A = (std::abs(numerator) > std::numeric_limits<double>::epsilon() or
                std::abs(denominator) > std::numeric_limits<double>::epsilon())
                   ? G * (numerator / denominator)
                   : 2.0 * G / (C + B);

    double P = (A * GAS_CONSTANT * T) / (pow(z, 2.0) * pow(FARADAY, 2.0));

    if (std::isnan(P)) {
        ArgErrLog("Failed to find permeability, check parameters");
    }

    return P;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurrent(double P, double V, int z, double T, double iconc, double oconc) {
    AssertLog(z != 0);
    AssertLog(T >= 0.0);
    AssertLog(iconc >= 0.0);
    AssertLog(oconc >= 0.0);

    double nuFoRT = static_cast<double>(z) * V * FARADAY / (GAS_CONSTANT * T);
    if (nuFoRT >= std::numeric_limits<double>::max_exponent10 * 2.30 ||
        nuFoRT <= std::numeric_limits<double>::min_exponent10 * 2.30) {
        throw std::runtime_error("Overflow encountered, nuFoRT: " + std::to_string(nuFoRT));
    }
    double eNuFoRT = std::exp(-nuFoRT);

    if (std::abs(nuFoRT) > std::numeric_limits<double>::epsilon()) {
        return P * z * FARADAY * nuFoRT * (iconc - oconc * eNuFoRT) / (1.0 - eNuFoRT);
    } else {
        // The limit of the indeterminate form is computed with L’Hôpital’s rule
        return P * z * FARADAY * (iconc - oconc);
    }
}

}  // namespace steps::math
