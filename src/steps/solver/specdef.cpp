/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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


/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

// STL headers.
#include <cassert>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/spec.hpp"
#include "steps/solver/specdef.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/types.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef::Specdef(Statedef * sd, uint idx, steps::model::Spec * s)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pSetupdone(false)
{
    AssertLog(pStatedef != nullptr);
    AssertLog(s != nullptr);
    pName = s->getID();

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef::~Specdef()
= default;

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::checkpoint(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::restore(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Specdef::name() const
{
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::setup()
{
  pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////


// END

