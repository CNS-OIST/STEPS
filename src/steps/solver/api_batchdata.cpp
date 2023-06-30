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

// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "api.hpp"
#include "compdef.hpp"
#include "patchdef.hpp"
#include "specdef.hpp"
#include "statedef.hpp"
// utils
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

std::vector<double> API::getBatchTetSpecCounts(const std::vector<index_t>& /* tets */,
                                               std::string const& /* s */) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getBatchTetSpecConcs(const std::vector<index_t>& /* tets */,
                                              std::string const& /* s */) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setBatchTetSpecConcs(const std::vector<index_t>& /*tets*/,
                               std::string const& /*s*/,
                               const std::vector<double>& /*concs*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getBatchTriSpecCounts(const std::vector<index_t>& /* tris */,
                                               std::string const& /* s */) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTetSpecCountsNP(const index_t* /* indices */,
                                  size_t /* input_size */,
                                  std::string const& /* s */,
                                  double* /* counts */,
                                  size_t /* output_size */) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTetSpecConcsNP(const index_t* /* indices */,
                                 size_t /* input_size */,
                                 std::string const& /* s */,
                                 double* /* counts */,
                                 size_t /* output_size */) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setBatchTetSpecConcsNP(const index_t* /*indices*/,
                                 size_t /*ntets*/,
                                 std::string const& /*s*/,
                                 const double* /*concs*/,
                                 size_t /*output_size*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTriSpecCountsNP(const index_t* /* indices */,
                                  size_t /* input_size */,
                                  std::string const& /* s */,
                                  double* /* counts */,
                                  size_t /* output_size */) const {
    NotImplErrLog("");
}

}  // namespace steps::solver
