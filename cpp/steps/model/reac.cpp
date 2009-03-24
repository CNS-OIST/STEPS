////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/error.hpp>
#include <steps/model/model.hpp>
#include <steps/model/volsys.hpp>
#include <steps/model/reac.hpp>
#include <steps/model/spec.hpp>

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

Reac::Reac(string const & id, Volsys * volsys, vector<Spec *> const & lhs,
           vector<Spec *> const & rhs, double kcst)
: pID(id)
, pModel(0)
, pVolsys(volsys)
, pLHS()
, pRHS()
, pOrder(0)
, pKcst(kcst)
{
    if (pVolsys == 0)
    {
        ostringstream os;
        os << "No volsys provided to Reac initializer function";
        throw steps::ArgErr(os.str());
    }
	if (pKcst < 0.0)
	{
		ostringstream os;
		os << "Reaction constant can't be negative";
		throw steps::ArgErr(os.str());
	}

	pModel = pVolsys->getModel();
	assert (pModel != 0);

	setLHS(lhs);
    setRHS(rhs);

    pVolsys->_handleReacAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Reac::~Reac(void)
{
    if (pVolsys == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Reac::_handleSelfDelete(void)
{
	pVolsys->_handleReacDel(this);
	pKcst = 0.0;
	pOrder = 0;
	pRHS.clear();
	pLHS.clear();
	pVolsys = 0;
	pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setID(string const & id)
{
    assert(pVolsys != 0);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVolsys->_handleReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setLHS(vector<Spec *> const & lhs)
{
    assert(pVolsys != 0);
    pLHS.clear();

    SpecPVecCI l_end = lhs.end();
    for (SpecPVecCI l = lhs.begin(); l != l_end; ++l)
    {
        assert ((*l)->getModel() == pModel);
        pLHS.push_back(*l);
    }
    pOrder = pLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setRHS(vector<Spec *> const & rhs)
{
    assert (pVolsys != 0);
    pRHS.clear();

    SpecPVecCI r_end = rhs.end();
    for (SpecPVecCI r = rhs.begin(); r != r_end; ++r)
    {
		assert ((*r)->getModel() == pModel);
        pRHS.push_back(*r);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setKcst(double kcst)
{
	assert(pVolsys != 0);
	if (kcst < 0.0)
	{
		ostringstream os;
		os << "Reaction constant can't be negative";
		throw steps::ArgErr(os.str());
	}
	pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

vector<Spec *> Reac::getAllSpecs(void) const
{
    SpecPVec specs = SpecPVec();
	bool first_occ = true;

	SpecPVec lhs = getLHS();
	SpecPVecCI l_end = lhs.end();
	for (SpecPVecCI l = lhs.begin(); l != l_end; ++l)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*l))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*l));
	}

	SpecPVec rhs = getRHS();
	SpecPVecCI r_end = rhs.end();
	for (SpecPVecCI r = rhs.begin(); r != r_end; ++r)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*r))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*r));
	}

	return specs;
}

////////////////////////////////////////////////////////////////////////////////

// END
