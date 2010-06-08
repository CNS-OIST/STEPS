////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
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

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

// STL headers.
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "geom.hpp"
#include "comp.hpp"
#include "patch.hpp"



////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
NAMESPACE_ALIAS(steps::wm, swm);

////////////////////////////////////////////////////////////////////////////////

swm::Comp::Comp(string const & id, swm::Geom * container, double vol)
: pID(id)
, pContainer(container)
, pVolsys()
, pVol(vol)
, pIPatches()
, pOPatches()
{
	if (pContainer == 0)
    {
        ostringstream os;
        os << "No container provided to Comp initializer function";
        throw steps::ArgErr(os.str());
    }

	if (pVol < 0.0)
	{
		ostringstream os;
		os << "Compartment volume can't be negative";
		throw steps::ArgErr(os.str());
	}
	pContainer->_handleCompAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

swm::Comp::~Comp(void)
{
	if (pContainer == 0) return;
	_handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::setID(string const & id)
{
	assert(pContainer != 0);
	if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer->_handleCompIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::setVol(double vol)
{
	assert(pContainer != 0);
	if (vol < 0.0)
	{
		ostringstream os;
		os << "Compartment volume can't be negative";
		throw steps::ArgErr(os.str());
	}
	pVol = vol;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::addVolsys(std::string const & id)
{
	// string identifier is only added to set if it is not already included
	pVolsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::delVolsys(std::string const & id)
{
	// string identifier is only removed from set if it is included
	pVolsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_addIPatch(swm::Patch * patch)
{
	assert (patch->getOComp() == this);
	// patch pointer is only added to set if it is not already included
	pIPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_delIPatch(swm::Patch * patch)
{
	assert (patch->getOComp() == this);
	pIPatches.erase(patch);
}
////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_addOPatch(swm::Patch * patch)
{
	assert (patch->getIComp() == this);
	// patch pointer is only added to set if it is not already included
	pOPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_delOPatch(swm::Patch * patch)
{
	assert (patch->getIComp() == this);
	pOPatches.erase(patch);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_handleSelfDelete(void)
{
	pContainer->_handleCompDel(this);
	pVol = 0.0;
	pVolsys.clear();
	pIPatches.clear();
	pOPatches.clear();
	pContainer = 0;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * swm::Comp::_getIPatch(uint lidx) const
{
	assert(lidx < pIPatches.size());
	std::set<steps::wm::Patch *>::const_iterator pit = pIPatches.begin();
	for (uint i=0; i < lidx; ++i) ++pit;
	return (*pit);
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * swm::Comp::_getOPatch(uint lidx) const
{
	assert(lidx < pOPatches.size());
	std::set<steps::wm::Patch *>::const_iterator pit = pOPatches.begin();
	for (uint i=0; i < lidx; ++i) ++pit;
	return (*pit);
}

////////////////////////////////////////////////////////////////////////////////

/// END
