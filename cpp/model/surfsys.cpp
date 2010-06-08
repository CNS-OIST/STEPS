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
#include <map>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "model.hpp"
#include "spec.hpp"
#include "surfsys.hpp"
#include "sreac.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

Surfsys::Surfsys(string const & id, Model * model)
: pID(id)
, pModel(model)
, pSReacs()
{
    if (pModel == 0)
    {
        ostringstream os;
        os << "No model provided to Surfsys initializer function";
        throw steps::ArgErr(os.str());
    }
    pModel->_handleSurfsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Surfsys::~Surfsys(void)
{
    if (pModel == 0) return;
	_handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::setID(string const & id)
{
    assert(pModel != 0);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleSurfsysIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSelfDelete(void)
{
	std::vector<steps::model::SReac *> allsreacs = getAllSReacs();
	SReacPVecCI sreac_end = allsreacs.end();
	for(SReacPVecCI sreac = allsreacs.begin(); sreac != sreac_end; ++sreac)
	{
		delete(*sreac);
	}

	pModel->_handleSurfsysDel(this);
	pSReacs.clear();
	pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

SReac * Surfsys::getSReac(string const & id) const
{
    SReacPMapCI sreac = pSReacs.find(id);
    if (sreac == pSReacs.end())
    {
        ostringstream os;
        os << "Model does not contain surface "
        "reaction with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(sreac->second != 0);
    return sreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delSReac(string const & id)
{
    SReac * sreac = getSReac(id);
	delete(sreac);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<SReac *> Surfsys::getAllSReacs(void) const
{
	SReacPVec sreacs = SReacPVec();
	SReacPMapCI sr_end = pSReacs.end();
	for (SReacPMapCI sr = pSReacs.begin(); sr != sr_end; ++sr)
	{
		sreacs.push_back(sr->second);
	}
	return sreacs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Surfsys::getAllSpecs(void) const
{
	SpecPVec specs = SpecPVec();
	bool first_occ = true;

	SReacPVec sreacs = getAllSReacs();
	SReacPVecCI sreac_end = sreacs.end();
	for (SReacPVecCI sreac = sreacs.begin(); sreac != sreac_end; ++sreac)
	{
		SpecPVec sr_specs = (*sreac)->getAllSpecs();
		SpecPVecCI sr_spec_end = sr_specs.end();
		for (SpecPVecCI sr_spec = sr_specs.begin();
		     sr_spec != sr_spec_end; ++sr_spec)
		{
			first_occ = true;
			SpecPVecCI allspecs_end = specs.end();
			for (SpecPVecCI allspecs = specs.begin();
				allspecs != allspecs_end; ++allspecs)
			{
				if ((*sr_spec) == (*allspecs))
				{
					first_occ = false;
					break;
				}
			}
			if (first_occ == true) specs.push_back(*sr_spec);
		}
	}

	return specs;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkSReacID(string const & id) const
{
    checkID(id);
    if (pSReacs.find(id) != pSReacs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacIDChange(string const & o, string const & n)
{
    SReacPMapCI sr_old = pSReacs.find(o);
    assert(sr_old != pSReacs.end());

    if(o==n) return;
    _checkSReacID(n);

    SReac * sr = sr_old->second;
    assert(sr != 0);
    pSReacs.erase(sr->getID());
    pSReacs.insert(SReacPMap::value_type(n,sr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacAdd(SReac * sreac)
{
    assert(sreac->getSurfsys() == this);
    _checkSReacID(sreac->getID());
    pSReacs.insert(SReacPMap::value_type(sreac->getID(), sreac));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacDel(SReac * sreac)
{
	assert (sreac->getSurfsys() == this);
    pSReacs.erase(sreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSpecDelete(Spec * spec)
{
	SReacPMapCI sreac_end = pSReacs.end();
	std::vector<std::string> sreacs_del = std::vector<std::string>();
    for (SReacPMapCI sreac = pSReacs.begin(); sreac != sreac_end; ++sreac)
    {
		SpecPVec specs = (sreac->second->getAllSpecs());
		SpecPVecCI sr_spec_end = specs.end();
		for (SpecPVecCI sr_spec = specs.begin();
			 sr_spec != sr_spec_end; ++sr_spec)
		{
			if ((*sr_spec) == spec)
			{
				sreacs_del.push_back(sreac->second->getID());
				break;
			}
		}
	}

	std::vector<std::string>::const_iterator sr_del_end = sreacs_del.end();
	for (std::vector<std::string>::const_iterator sr_del = sreacs_del.begin();
		 sr_del != sr_del_end; ++sr_del)
	{
		delSReac(*sr_del);
	}
}

////////////////////////////////////////////////////////////////////////////////

SReac * Surfsys::_getSReac(uint lidx) const
{
	assert (lidx < pSReacs.size());
	std::map<std::string, SReac *>::const_iterator sr_it = pSReacs.begin();
	for (uint i=0; i< lidx; ++i) ++sr_it;
	return sr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

// END
