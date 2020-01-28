/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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

std::vector<double> API::getROITetCounts(const std::string& /*ROI_id*/, std::string const & /*s*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITriCounts(const std::string& /*ROI_id*/, std::string const & /*s*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITetCountsNP(const std::string& /*ROI_id*/, std::string const & /*s*/, double* /*counts*/, int /*output_size*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITriCountsNP(const std::string& /*ROI_id*/, std::string const & /*s*/, double* /*counts*/, int /*output_size*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIVol(const std::string& /*ROI_id*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIArea(const std::string& /*ROI_id*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROICount(const std::string& /*ROI_id*/, std::string const & /*s*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROICount(const std::string& /*ROI_id*/, std::string const & /*s*/, double /*count*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIAmount(const std::string& /*ROI_id*/, std::string const & /*s*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIConc(const std::string& /*ROI_id*/, std::string const & /*s*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIConc(const std::string& /*ROI_id*/, std::string const & /*s*/, double /*conc*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIClamped(const std::string& /*ROI_id*/, std::string const & /*s*/, bool /*b*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIReacK(const std::string& /*ROI_id*/, std::string const & /*r*/, double /*kf*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROISReacK(const std::string& /*ROI_id*/, std::string const & /*sr*/, double /*kf*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIDiffD(const std::string& /*ROI_id*/, std::string const & /*d*/, double /*dk*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIReacActive(const std::string& /*ROI_id*/, std::string const & /*r*/, bool /*a*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROISReacActive(const std::string& /*ROI_id*/, std::string const & /*sr*/, bool /*a*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIDiffActive(const std::string& /*ROI_id*/, std::string const & /*d*/, bool /*act*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIVDepSReacActive(const std::string& /*ROI_id*/, std::string const & /*vsr*/, bool /*a*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getROIReacExtent(const std::string& /*ROI_id*/, std::string const & /*r*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROIReacExtent(const std::string& /*ROI_id*/, std::string const & /*r*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getROISReacExtent(const std::string& /*ROI_id*/, std::string const & /*sr*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROISReacExtent(const std::string& /*ROI_id*/, std::string const & /*sr*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getROIDiffExtent(const std::string& /*ROI_id*/, std::string const & /*d*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROIDiffExtent(const std::string& /*ROI_id*/, std::string const & /*d*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////


// END

