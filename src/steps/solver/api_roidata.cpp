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

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITetCounts(std::string ROI_id, std::string const & s) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITriCounts(std::string ROI_id, std::string const & s) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITetCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITriCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIVol(std::string ROI_id) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIArea(std::string ROI_id) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getROICount(std::string ROI_id, std::string const & s) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROICount(std::string ROI_id, std::string const & s, double count)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIAmount(std::string ROI_id, std::string const & s) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIConc(std::string ROI_id, std::string const & s) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIConc(std::string ROI_id, std::string const & s, double conc)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIClamped(std::string ROI_id, std::string const & s, bool b)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIReacK(std::string ROI_id, std::string const & r, double kf)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROISReacK(std::string ROI_id, std::string const & sr, double kf)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIDiffD(std::string ROI_id, std::string const & d, double dk)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIReacActive(std::string ROI_id, std::string const & r, bool a)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROISReacActive(std::string ROI_id, std::string const & sr, bool a)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIDiffActive(std::string ROI_id, std::string const & d, bool act)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIVDepSReacActive(std::string ROI_id, std::string const & vsr, bool a)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getROIReacExtent(std::string ROI_id, std::string const & r) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROIReacExtent(std::string ROI_id, std::string const & r)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getROISReacExtent(std::string ROI_id, std::string const & sr) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROISReacExtent(std::string ROI_id, std::string const & sr)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getROIDiffExtent(std::string ROI_id, std::string const & d) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::resetROIDiffExtent(std::string ROI_id, std::string const & d)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////


// END

