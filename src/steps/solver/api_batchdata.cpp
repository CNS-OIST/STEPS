/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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


// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/api.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/specdef.hpp"
#include "steps/solver/statedef.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getBatchTetCounts(std::vector<uint> const & tets, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getBatchTriCounts(std::vector<uint> const & tris, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END

