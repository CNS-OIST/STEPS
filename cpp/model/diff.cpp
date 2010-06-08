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
#include "model.hpp"
#include "volsys.hpp"
#include "diff.hpp"
#include "spec.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(string const & id, Volsys * volsys, Spec * lig, double dcst)
: pID(id)
, pModel(0)
, pVolsys(volsys)
, pLig(lig)
, pDcst(dcst)
{
    if (pVolsys == 0)
    {
        ostringstream os;
        os << "No volsys provided to Diff initializer function.";
        throw steps::ArgErr(os.str());
    }
	if(pDcst < 0.0)
	{
		ostringstream os;
		os << "Diffusion constant can't be negative";
		throw steps::ArgErr(os.str());
	}
    pModel = pVolsys->getModel();
    assert (pModel != 0);

    pVolsys->_handleDiffAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::~Diff(void)
{
    if (pVolsys == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Diff::_handleSelfDelete(void)
{
	pVolsys->_handleDiffDel(this);
	pDcst = 0.0;
	pLig = 0;
	pVolsys = 0;
	pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setID(string const & id)
{
    assert(pVolsys != 0);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVolsys->_handleDiffIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDcst(double dcst)
{
	assert(pVolsys != 0);
	if(dcst < 0.0)
	{
		ostringstream os;
		os << "Diffusion constant can't be negative";
		throw steps::ArgErr(os.str());
	}
	pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setLig(Spec * lig)
{
    assert(lig !=0);
    pLig = lig;
}

////////////////////////////////////////////////////////////////////////////////

vector<Spec *> Diff::getAllSpecs(void) const
{
    SpecPVec specs = SpecPVec();
	specs.push_back(pLig);
	return specs;
}

////////////////////////////////////////////////////////////////////////////////

// END
