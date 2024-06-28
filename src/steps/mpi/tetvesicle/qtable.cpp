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


// BOOST
#include <boost/math/quadrature/gauss_kronrod.hpp>

// STEPS headers.
#include "math/constants.hpp"  //PI

namespace steps::mpi::tetvesicle {

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

void Qtable::setup() noexcept {
    // Do the integration to pi to get the normalisation constant
    double result, error;
    std::function<double(double)> f = std::bind(&Qtable::Q, this, std::placeholders::_1);

    result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
        f, 0, math::PI, 5, 1e-7, &error);

    double norm_constant = result;


    // Use a large table to set up the cdf table
    const unsigned int nvals = 10001;
    double x[nvals];
    double y_cdf[nvals];

    for (unsigned int i = 0; i < nvals; ++i) {
        double xval = i * (math::PI / (nvals - 1));
        x[i] = xval;

        result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
            f, 0, xval, 5, 1e-7, &error);

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
        pQ_values[i] = Q(xinterp);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Qtable::reinit(unsigned int size, double tau) noexcept {
    pTablesize = size;
    pTau = tau;

    pX_interp.resize(pTablesize + 1);
    pQ_values.resize(pTablesize + 1);

    setup();
}

////////////////////////////////////////////////////////////////////////////////

double Qtable::getPhi() const noexcept {
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
            if (y < Q(xval)) {
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

////////////////////////////////////////////////////////////////////////////////

double Qtable::Q(double theta) const noexcept {
    return (1.0 / pTau) * std::sqrt(theta * std::sin(theta)) *
           std::exp(-(std::pow(theta, 2)) / (2 * pTau));
}

}  // namespace steps::mpi::tetvesicle
