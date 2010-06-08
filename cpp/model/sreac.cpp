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
#include <iostream>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "model.hpp"
#include "surfsys.hpp"
#include "sreac.hpp"
#include "spec.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

SReac::SReac(string const & id, Surfsys * surfsys,
             vector<Spec *> const & olhs, vector<Spec *> const & ilhs,
             vector<Spec *> const & slhs,
             vector<Spec *> const & irhs, vector<Spec *> const & srhs,
             vector<Spec *> const & orhs, double kcst)
: pID(id)
, pModel(0)
, pSurfsys(surfsys)
, pOuter(true)
, pOLHS()
, pILHS()
, pSLHS()
, pIRHS()
, pSRHS()
, pORHS()
, pOrder(0)
, pKcst(kcst)
{
    if (pSurfsys == 0)
    {
        ostringstream os;
        os << "No surfsys provided to SReac initializer function";
        throw steps::ArgErr(os.str());
    }
	if (pKcst < 0.0)
	{
		ostringstream os;
		os << "Surface reaction constant can't be negative";
		throw steps::ArgErr(os.str());
	}

    // Can't have species on the lhs in the inner and outer compartment
    if (olhs.size() != 0 && ilhs.size() != 0)
    {
    	ostringstream os;
    	os << "Volume lhs species must belong to either inner or outer ";
    	os << "compartment, not both.";
    	throw steps::ArgErr(os.str());
    }

	pModel = pSurfsys->getModel();
	assert (pModel != 0);

    if (olhs.size() > 0) setOLHS(olhs);
    if (ilhs.size() > 0) setILHS(ilhs);
    setSLHS(slhs);
    setIRHS(irhs);
    setSRHS(srhs);
    setORHS(orhs);

    pSurfsys->_handleSReacAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

SReac::~SReac(void)
{
    if (pSurfsys == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::_handleSelfDelete(void)
{
	pSurfsys->_handleSReacDel(this);
	pKcst = 0.0;
	pOrder = 0;
	pORHS.clear();
	pSRHS.clear();
	pIRHS.clear();
	pSLHS.clear();
	pILHS.clear();
	pOLHS.clear();
	pSurfsys = 0;
	pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setID(string const & id)
{
    assert(pSurfsys != 0);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys->_handleSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setOLHS(vector<Spec *> const & olhs)
{
    assert (pSurfsys != 0);

	if (pILHS.size() != 0)
	{
    	ostringstream os;
    	os << "WARNING: Removing inner compartment species from lhs stoichiometry";
	}
    pILHS.clear();
	pOLHS.clear();
    SpecPVecCI ol_end = olhs.end();
    for (SpecPVecCI ol = olhs.begin(); ol != ol_end; ++ol)
    {
		assert ((*ol)->getModel() == pModel);
        pOLHS.push_back(*ol);
    }
    pOuter = true;
    pOrder = pOLHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setILHS(vector<Spec *> const & ilhs)
{
    assert (pSurfsys != 0);

	if (pOLHS.size() != 0)
	{
    	ostringstream os;
    	os << "WARNING: Removing outer compartment species from lhs stoichiometry";
	}
    pOLHS.clear();
	pILHS.clear();
    SpecPVecCI il_end = ilhs.end();
    for (SpecPVecCI il = ilhs.begin(); il != il_end; ++il)
    {
		assert ((*il)->getModel() == pModel);
        pILHS.push_back(*il);
    }
    pOuter = false;
    pOrder = pILHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setSLHS(vector<Spec *> const & slhs)
{
    assert (pSurfsys != 0);
    pSLHS.clear();
    SpecPVecCI sl_end = slhs.end();
    for (SpecPVecCI sl = slhs.begin(); sl != sl_end; ++sl)
    {
		assert ((*sl)->getModel() == pModel);
        pSLHS.push_back(*sl);
    }

    if (pOuter) { pOrder = pOLHS.size() + pSLHS.size(); }
    else { pOrder = pILHS.size() + pSLHS.size(); }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setIRHS(vector<Spec *> const & irhs)
{
    assert (pSurfsys != 0);
    pIRHS.clear();
    SpecPVecCI ir_end = irhs.end();
    for (SpecPVecCI ir = irhs.begin(); ir != ir_end; ++ir)
    {
		assert ((*ir)->getModel() == pModel);
        pIRHS.push_back(*ir);
    }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setSRHS(vector<Spec *> const & srhs)
{
    assert (pSurfsys != 0);
    pSRHS.clear();
    SpecPVecCI sr_end = srhs.end();
    for (SpecPVecCI sr = srhs.begin(); sr != sr_end; ++sr)
    {
        assert ((*sr)->getModel() == pModel);
        pSRHS.push_back(*sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setORHS(vector<Spec *> const & orhs)
{
    assert (pSurfsys != 0);
    pORHS.clear();
    SpecPVecCI or_end = orhs.end();
    for (SpecPVecCI ors = orhs.begin(); ors != or_end; ++ors)
    {
		assert ((*ors)->getModel() == pModel);
        pORHS.push_back(*ors);
    }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setKcst(double kcst)
{
	assert (pSurfsys != 0);
	if(kcst < 0.0)
	{
		ostringstream os;
		os << "Surface reaction constant can't be negative";
		throw steps::ArgErr(os.str());
	}
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

vector<Spec *> SReac::getAllSpecs(void) const
{
    SpecPVec specs = SpecPVec();
	bool first_occ = true;
	assert(pOLHS.size() == 0 || pILHS.size() == 0);

	SpecPVec olhs = getOLHS();
	SpecPVecCI ol_end = olhs.end();
	for (SpecPVecCI ol = olhs.begin(); ol != ol_end; ++ol)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*ol))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*ol));
	}

	SpecPVec ilhs = getILHS();
	SpecPVecCI il_end = ilhs.end();
	for (SpecPVecCI il = ilhs.begin(); il != il_end; ++il)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*il))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*il));
	}

	SpecPVec slhs = getSLHS();
	SpecPVecCI sl_end = slhs.end();
	for (SpecPVecCI sl = slhs.begin(); sl != sl_end; ++sl)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*sl))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*sl));
	}

	SpecPVec irhs = getIRHS();
	SpecPVecCI ir_end = irhs.end();
	for (SpecPVecCI ir = irhs.begin(); ir != ir_end; ++ir)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*ir))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*ir));
	}

	SpecPVec srhs = getSRHS();
	SpecPVecCI sr_end = srhs.end();
	for (SpecPVecCI sr = srhs.begin(); sr != sr_end; ++sr)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*sr))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*sr));
	}

	SpecPVec orhs = getORHS();
	SpecPVecCI ors_end = orhs.end();
	for (SpecPVecCI ors = orhs.begin(); ors != ors_end; ++ors)
	{
		first_occ = true;
		SpecPVecCI s_end = specs.end();
		for (SpecPVecCI s = specs.begin(); s != s_end; ++s)
		{
			if ((*s) == (*ors))
			{
				first_occ = false;
				break;
			}
		}
		if (first_occ == true) specs.push_back((*ors));
	}

	return specs;
}

////////////////////////////////////////////////////////////////////////////////

// END
