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

#include "mpi/tetvesicle/qtable.hpp"

// Standard library & STL headers.
#include <chrono>
#include <cmath>
#include <random>

// logging

// GSL
#include <gsl/gsl_integration.h>

// STEPS headers.
#include "math/constants.hpp"  //PI

namespace steps::mpi::tetvesicle {

// For some reason this would only compile with the function definition outside
// of the class
double Q(double theta, void* params) {
    double tau = *(double*) params;
    return (1.0 / tau) * sqrt(theta * sin(theta)) * exp(-(pow(theta, 2)) / (2 * tau));
}

////////////////////////////////////////////////////////////////////////////////

Qtable::Qtable(unsigned int size, double tau, const rng::RNGptr& r)
    : pTablesize(size)
    , pTau(tau)
    , pX_interp(pTablesize + 1)
    , pQ_values(pTablesize + 1)
    , rng(r) {
    setup();
}

////////////////////////////////////////////////////////////////////////////////

void Qtable::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve. Qtables can only be created with API calls (setVesicleSpecDiffD,
    // setVesicleSurfaceLinkSpecSDiffD) so nothing to do here. Solver has to
    // take care of it.
}

////////////////////////////////////////////////////////////////////////////////

void Qtable::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Qtable::setup() {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    // Create the GSL integration function
    gsl_function F;
    F.function = &Q;
    F.params = &pTau;

    // Do the integration to pi to get the normalisation constant
    double result, error;
    gsl_integration_qags(&F, 0, math::PI, 0, 1e-7, 1000, w, &result, &error);
    double norm_constant = result;

    // Use a large table to set up the cdf table
    const unsigned int nvals = 10001;
    double x[nvals];      // to avoid compile warnings with double x[nvals]
    double y_cdf[nvals];  // to avoid compile warnings with double y_cdf[nvals]

    for (unsigned int i = 0; i < nvals; ++i) {
        double xval = i * (math::PI / (nvals - 1));
        x[i] = xval;

        gsl_integration_qags(&F, 0, xval, 0, 1e-7, 1000, w, &result, &error);

        y_cdf[i] = result / norm_constant;
    }

    pX_interp[pTablesize] = math::PI;
    pQ_values[pTablesize] = 0;

    // setup the x (theta) interpolation table at constant steps of y (cdf)
    unsigned int idx = 0;
    for (unsigned int i = 0; i < pTablesize; ++i) {
        double y_val = i / static_cast<double>(pTablesize);

        while (y_cdf[idx + 1] < y_val) {
            idx += 1;
        }

        double xinterp = x[idx] + ((x[idx + 1] - x[idx]) / (y_cdf[idx + 1] - y_cdf[idx])) *
                                      (y_val - y_cdf[idx]);

        pX_interp[i] = xinterp;
        pQ_values[i] = Q(xinterp, &pTau);
    }

    gsl_integration_workspace_free(w);
}

////////////////////////////////////////////////////////////////////////////////

void Qtable::reinit(unsigned int size, double tau) {
    pTablesize = size;
    pTau = tau;

    pX_interp.resize(pTablesize + 1);
    pQ_values.resize(pTablesize + 1);

    setup();
}

////////////////////////////////////////////////////////////////////////////////

double Qtable::getPhi() {
    constexpr unsigned int maxRetries = 10000;
    constexpr double threshold = 0.1;

    double rnum = pTablesize * rng->getUnfII();

    // Find xvalue:
    int idx = static_cast<int>(rnum);
    double rfrac = (rnum - idx);
    double xval = (1 - rfrac) * pX_interp[idx] + rfrac * pX_interp[idx + 1];

    // Use rejection sampling if necessary
    if (pX_interp[idx + 1] - pX_interp[idx] > threshold) {
        unsigned int retries = 0;
        double maxQ = std::max(pQ_values[idx], pQ_values[idx + 1]);
        while (retries < maxRetries) {
            double y = maxQ * rng->getUnfII();
            if (y < Q(xval, &pTau)) {
                break;
            } else {
                ++retries;
                rfrac = rng->getUnfII();
                xval = (1 - rfrac) * pX_interp[idx] + rfrac * pX_interp[idx + 1];
            }
        }
        if (retries == maxRetries) {
            // If all samples were rejected, take the lowest value
            xval = pX_interp[idx];
        }
    }

    return xval;
}

}  // namespace steps::mpi::tetvesicle
