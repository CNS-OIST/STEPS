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
#include <string>
#include <sstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/specdef.hpp"
// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITetCounts(std::string ROI_id, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITriCounts(std::string ROI_id, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITetCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITriCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIVol(std::string ROI_id) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIArea(std::string ROI_id) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROICount(std::string ROI_id, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROICount(std::string ROI_id, std::string const & s, double count)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIAmount(std::string ROI_id, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIConc(std::string ROI_id, std::string const & s) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIConc(std::string ROI_id, std::string const & s, double conc)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIClamped(std::string ROI_id, std::string const & s, bool b)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIReacK(std::string ROI_id, std::string const & r, double kf)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROISReacK(std::string ROI_id, std::string const & sr, double kf)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIDiffD(std::string ROI_id, std::string const & d, double dk)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIReacActive(std::string ROI_id, std::string const & r, bool a)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROISReacActive(std::string ROI_id, std::string const & sr, bool a)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIDiffActive(std::string ROI_id, std::string const & d, bool act)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIVDepSReacActive(std::string ROI_id, std::string const & vsr, bool a)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

uint API::getROIReacExtent(std::string ROI_id, std::string const & r) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROIReacExtent(std::string ROI_id, std::string const & r)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

uint API::getROISReacExtent(std::string ROI_id, std::string const & sr) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROISReacExtent(std::string ROI_id, std::string const & sr)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

uint API::getROIDiffExtent(std::string ROI_id, std::string const & d) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROIDiffExtent(std::string ROI_id, std::string const & d)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////


// END

