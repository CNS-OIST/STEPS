////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2011ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////


// STL headers.
#include <string>
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "api.hpp"
#include "statedef.hpp"

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

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

void API::_setDiffBoundaryDiffusionActive(uint dbidx, uint sidx, bool act)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getDiffBoundaryDiffusionActive(uint dbidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END
