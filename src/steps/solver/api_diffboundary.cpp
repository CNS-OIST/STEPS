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



// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "api.hpp"
#include "statedef.hpp"
// util
#include "util/error.hpp"
// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

void API::setDiffBoundaryDiffusionActive(string const & db, string const & s, bool act)
{
    uint dbidx = pStatedef->getDiffBoundaryIdx(db);
    uint sidx = pStatedef->getSpecIdx(s);

    return _setDiffBoundaryDiffusionActive(dbidx, sidx, act);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getDiffBoundaryDiffusionActive(string const & db, string const & s) const
{
    uint dbidx = pStatedef->getDiffBoundaryIdx(db);
    uint sidx = pStatedef->getSpecIdx(s);

    return _getDiffBoundaryDiffusionActive(dbidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setDiffBoundaryDcst(std::string const & db, std::string const & s, double dcst, std::string const & direction_comp)
{
    uint dbidx = pStatedef->getDiffBoundaryIdx(db);
    uint sidx = pStatedef->getSpecIdx(s);
    if (direction_comp.empty()) {
        _setDiffBoundaryDcst(dbidx, sidx, dcst);
    }
    else {
        uint cidx = pStatedef->getCompIdx(direction_comp);
        _setDiffBoundaryDcst(dbidx, sidx, dcst, cidx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::_setDiffBoundaryDiffusionActive(uint /*dbidx*/, uint /*sidx*/, bool /*act*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getDiffBoundaryDiffusionActive(uint /*dbidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setDiffBoundaryDcst(uint /*dbidx*/, uint /*sidx*/, double /*dcst*/, uint /*direction_comp*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END
